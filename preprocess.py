# coding=UTF-8


'''
import module
'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
import shlex, subprocess
from subprocess import Popen,PIPE
import time
from matplotlib import cm
import seaborn as sns
from datetime import datetime
import argparse
import subprocess
import pyBigWig
import time
import os, sys


def get_parser():
    """
    Create a parser and add arguments
    """
    parser = argparse.ArgumentParser(usage="Usage: preprocess.py {-s <samplelist> | -n <samplename> -f <cgmap> -e <expressionfile>} -g <gtf> [options]") 
    group2 = parser.add_argument_group('For samplelist')
    group2.add_argument("-s","--samplelist",type=str,default="False",help="input sample description file")
    group1 = parser.add_argument_group('For individual data')
    group1.add_argument("-n","--samplename",type=str,help="Setting the name of the set of data. This determines the names of output files for downstream analyses")
    group1.add_argument("-f","--cgmapfilename",type=str,help="input CGmap file")
    group1.add_argument("-e","--expressionfile",type=str,help="input gene expression file")
    group3 = parser.add_argument_group('General options')
    group3.add_argument("-g","--gtffile",type=str,help="input gene annotation file")
    group3.add_argument("-c","--cutoff",type=int,default=4, help="minimum cytosines that are covered by reads")
    return parser



def cgmap2bwbz(name,cgfile,cutoff):
    """
    cut the depth of cgmap, and save 
    """
    if cgfile[-3:]==".gz":
    	cgmap=pd.read_csv(cgfile,sep='\t',header=None,dtype={0:str},compression='gzip') 
    elif cgfile[-6:]==".CGmap":
    	cgmap=pd.read_csv(cgfile,sep="\t",header=None,dtype={0:str})
    else:
    	try:
    		cgmap=pd.read_csv(cgfile,sep="\t",header=None,dtype={0:str})
    	except:
    		print "wrong format of the input cgmap"
    		##若前面已經錯了，迴圈應該也不用跑了吧>>>怎麼做
    cgmap_cutoff=cgmap[cgmap[7]>cutoff]
    if len(cgmap_cutoff)>0:
    	print "finish cutoff"
    else:
    	print "no file"
    if cgfile[-3:]==".gz":
        cgmap_cutoff.to_csv("%s.cgmap.gz"%(name),header=False,sep="\t",index=False,compression='gzip') ##用成gzip比較久
        print "save file %s.cgmap.gz"%(name)
        bzfile=name+".bz"
        subprocess.call("cgmaptools convert cgmap2cgbz -c %s -b %s"%(cgfile,bzfile),shell=True)
        print "save %s"%(bzfile)
    elif cgfile[-6:]==".CGmap":
        cgmap_cutoff.to_csv("cutoff%d_%s.gz"%(cutoff,cgfile),header=False,sep="\t",index=False,compression='gzip')
        print "save file cutoff%d_%s.gz"%(cutoff,cgfile)
        bzfile=name+".bz"
        subprocess.call("cgmaptools convert cgmap2cgbz -c %s -b %s"%(cgfile,bzfile),shell=True)
        print "save %s"%(bzfile)
    cgmapc=cgmap_cutoff.copy()
    cgoutfile=name+"_CG.bw"
    chgoutfile=name+"_CHG.bw"
    chhoutfile=name+"_CHH.bw"
    ## the boundary of bw in each chr
    bwheader=[]
    for chr in cgmapc.iloc[:,0].unique():
    	k=cgmapc[cgmapc.iloc[:,0]==chr]
    	maxPos=k.iloc[:,2].max()#position max in the chr
    	bwheader.append((chr,maxPos+100000)) ##to make no default when run bw, no bound problem
    ## open bw file to write
    bw_CG = pyBigWig.open(cgoutfile,"w")
    bw_CG.addHeader(bwheader)
    bw_CHG = pyBigWig.open(chgoutfile,"w")
    bw_CHG.addHeader(bwheader)
    bw_CHH = pyBigWig.open(chhoutfile,"w")
    bw_CHH.addHeader(bwheader)
    ## write in bw info
    for chromosome in cgmapc.iloc[:,0].unique():
    	k=cgmapc[cgmapc.iloc[:,0]==chromosome]
    	CG = k[(k.iloc[:,3] == 'CG') ]#& (k.ix[:,5] != "-")
    	CHG = k[(k.iloc[:,3] == 'CHG')]#& (k.ix[:,5] != "-")
    	CHH = k[(k.iloc[:,3] == 'CHH')]#& (k.ix[:,5] != "-")
    	bw_CG.addEntries(chromosome,[x-1 for x in CG.ix[:,2].tolist()],values=[float(f) for f in CG.ix[:,5].tolist()], span=1)
    	bw_CHG.addEntries(chromosome,[x-1 for x in CHG.ix[:,2].tolist()],values=[float(f) for f in CHG.ix[:,5].tolist()], span=1)
    	bw_CHH.addEntries(chromosome,[x-1 for x in CHH.ix[:,2].tolist()],values=[float(f) for f in CHH.ix[:,5].tolist()], span=1)
    ##　close and finish the file
    bw_CG.close()
    bw_CHG.close()
    bw_CHH.close()
    bw_CG = pyBigWig.open(cgoutfile)
    bw_CHG = pyBigWig.open(chgoutfile)
    bw_CHH = pyBigWig.open(chhoutfile)
    bwall={}
    bwall["CG"]=bw_CG
    bwall["CHG"]=bw_CHG
    bwall["CHH"]=bw_CHH
    print "finish create bw file"
    return cgmap_cutoff,bwall,bzfile



def gtfToGenePred(gtf):  ##if add the sample name
    genePred="genePred"
    subprocess.call("gtfToGenePred %s genePred -genePredExt"%(gtf),shell=True)
    return genePred



def AllinfoAnnotation(genePred):
    subprocess.call("cat %s|awk \'OFS=\"\t\"{print $12,$1,$2,$3,$4,$5,$5-$4,$8,$9,$10}\'>all_info_annottion"%(genePred),shell=True)##$12genename>>from extend genePred
    df_a=pd.read_csv("all_info_annottion",sep='\t',header=None)
    big_bed = pd.DataFrame()
    gene_list=df_a[0].unique()
    for i in range(len(gene_list)):
        gene=df_a[df_a[0]==gene_list[i]] ## choose all data for that gene
        if len(gene)>1:
            ##can not change, choose the longest gene
            gene=gene.loc[gene[6]==max(gene[6])]
        if len(gene)>1:
            ##choose the more exons one
            gene=gene.loc[gene[7]==max(gene[7])]
        onlygene=gene.drop_duplicates(subset=[0], keep='first', inplace=False) ## if gene length & exon# is the same
        big_bed=big_bed.append(onlygene) ##.iloc should add>>  seperate from iloc first
    big_bed.to_csv("big_bed",sep="\t",index=False,header=False)
    print("save the information of gene annotation: big_bed")
    #big_bed=pd.read_csv("OTU5_%s_big_bed"%(args.date),sep="\t",header=None)
    print("Total gene number is %d"%(len(big_bed)))
    return big_bed


def Create_gb_bed(big_bed):
    '''
    gene body bed
    '''
    gb_bed=big_bed[[2,4,5,0,3]]  ##
    fill=np.full((len(gb_bed),1),".") ## np.full() 用於填充很方便, create a column
    gb_bed.insert(4,".",fill)## flle the fifth column
    #gbfile="OTU5_%s_gb_bed"%(args.date)
    gb_bed.to_csv("gb.bed",sep="\t",index=False,header=False)
    print("save Gene_Body_Bed")
    return gb_bed



def Create_pmt_bed():
    '''
    PMT bed
    '''
    subprocess.call("cat %s|awk \'OFS=\"\t\" {if  ($6==\"+\") print $1,$2-2000,$2,$4,$5,$6 ;else print $1,$3,$3+2000,$4,$5,$6}\'>pmt.bed"%("gb.bed"),shell=True)
    print("save Promoter_Bed")
    pmt_bed=pd.read_csv("pmt.bed",sep="\t",header=None,dtype={0:str})
    ##cat gb_bed|awk 'OFS="\t" {if  ($6=="+") print $1,$2-2000,$2,$4,$5,$6 ;else print $1,$3,$3+2000,$4,$5,$6}'>pmt_bed
    ##cat rpkm.txt|awk 'OFS="\t"{print $1,$2,$7} '|head   ## this is for rpkm table
    return pmt_bed

 
def CreateExonIntronBed(big_bed):
    exon_bed=pd.DataFrame()
    intron_bed=pd.DataFrame()
    exon_list = []
    intron_list = []
    exsta=list(big_bed[8])
    exend=list(big_bed[9])
    for i in range(len(big_bed)):
        exsta[i]=exsta[i].strip().split(",")
        exsta[i].pop()
        exend[i]=exend[i].strip().split(",")
        exend[i].pop()
        for j in range(len(exsta[i])): ##這兩行也許可以合併
            ##combine exon &intron bed
            exon_row=[big_bed.iloc[i,2],int(exsta[i][j]),int(exend[i][j]),big_bed.iloc[i,0],big_bed.iloc[i,3]] ###chr,sta,end,gene,dir
            #print exon_row
            exon_list.append(exon_row)
        for j in range(len(exsta[i])-1):
            intron_row=[big_bed.iloc[i,2],int(exend[i][j]),int(exsta[i][j+1]),big_bed.iloc[i,0],big_bed.iloc[i,3]]
            intron_list.append(intron_row)
    exon_bed = exon_bed.append(exon_list)
    intron_bed = intron_bed.append(intron_list)
    fill_exon=np.full((len(exon_bed),1),".") 
    exon_bed.insert(4,".",fill_exon)## flle the fifth column
    fill_intron=np.full((len(intron_bed),1),".") ## np.full() , create a column
    intron_bed.insert(4,".",fill_intron)## flle the fifth column
    exon_bed.to_csv("exon.bed",sep="\t",index=False,header=False)
    print("save exon bed")
    intron_bed.to_csv("intron.bed",sep="\t",index=False,header=False)
    print("save intron bed")
    return exon_bed,intron_bed



###average meth

def AverageMethGeneBW(bwall,chr,left,right):
    Meth_value_mean={}
    for i in ["CG","CHG","CHH"]:
        bw=bwall["%s"%(i)]
        Meth_mean=np.nanmean(bw.values(chr,left-1,right))
        Meth_value_mean["%s"%(i)]=Meth_mean
    return Meth_value_mean



def RunMethGeneBw(bwall,read_bed):
    df1 = pd.DataFrame()
    df1_list=[]
    for j in range(len(read_bed)):
        try:
            Meth_value_mean=AverageMethGeneBW(bwall,read_bed.iloc[j,0],read_bed.iloc[j,1],read_bed.iloc[j,2])
            s1 = list((read_bed.iloc[j,3],Meth_value_mean["CG"],Meth_value_mean["CHG"],Meth_value_mean["CHH"]))
            df1_list.append(s1)
            print j, s1
        except:
            print j  ###there are Pt
            pass
    df1 = df1.append(df1_list, ignore_index=True)
    df1.columns=['gene_ID','Meth_value_mean_CG',"Meth_value_mean_CHG","Meth_value_mean_CHH"]
    return df1



def func_Meth_siteC_siteCT(Meth_file,chr,left,right,cutoff):   ##calculate the number of C and the sum of all methylation level in the exon
    command_line ="cgmaptools fetch cgbz -b %s -C %s -L %d -R %d"%(Meth_file,chr,left,right)
    args = shlex.split(command_line)       
    a=subprocess.check_output(args, stdin=None, stderr=None, shell=False, universal_newlines=False)#print(a)#也是看起來是dataframe
    #b=a.replace("\t",",")#好像不用
    c=a.splitlines(True)#>>> 將\t的區別開來，分成很多個string(但還有\n)
    #d=c.split()這樣不行ㄚㄚㄚ因為大的是list
    #c[0].split()就可以
    for i in range(len(c)):
        c[i]=c[i].split()
    cgmap_select=pd.DataFrame(c,columns=["chr","1Letter","site","3letter","2letter","Meth_level","numC","numC_T"])
    cgmap_select=cgmap_select.convert_objects(convert_numeric=True)
    cgmap_select['chr'].astype(str)
    Meth_site_cut=cgmap_select[cgmap_select.numC_T>cutoff] ##cutoff first
    Meth_siteC_siteCT={}
    for i in ["CG","CHG","CHH"]:
        Meth_site_context=Meth_site_cut[Meth_site_cut["3letter"]==i]
        sumC=Meth_site_context["numC"].sum()
        sumCT=Meth_site_context["numC_T"].sum()
        Meth_siteC_siteCT["%s"%(i)]=[sumC,sumCT]
    return Meth_siteC_siteCT    ##可以Meth_siteC_siteCT["CG"][0]這樣取值


def RunMethExon(cgmap,extron_bed):
    ID_sumC_CT=pd.DataFrame() ## create all value for all exons in the bed. file 
    row_list=[]
    for i in range(len(extron_bed)):
        Meth_siteC_siteCT=func_Meth_siteC_siteCT(cgmap,extron_bed.iloc[i,0],extron_bed.iloc[i,1],extron_bed.iloc[i,2],4)
        row_value=[extron_bed.iloc[i,3],Meth_siteC_siteCT["CG"][0],Meth_siteC_siteCT["CG"][1],Meth_siteC_siteCT["CHG"][0],Meth_siteC_siteCT["CHG"][1],Meth_siteC_siteCT["CHH"][0],Meth_siteC_siteCT["CHH"][1]]
        row_list.append(row_value)
    ID_sumC_CT = ID_sumC_CT.append(row_list)
    ID_sumC_CT.columns=['gene_ID','siteC_CG',"siteCT_CG","siteC_CHG","siteCT_CHG","siteC_CHH","siteCT_CHH"]
    AverageMethExon = pd.DataFrame() ## calculate average methylation level of each gene
    mean_list = []
    for j in ID_sumC_CT["gene_ID"].unique():
        exons_gene=ID_sumC_CT.loc[ID_sumC_CT["gene_ID"]==j] ## choose all exons for that gene
        Meth_value_mean_CG=float(exons_gene["siteC_CG"].sum())/exons_gene["siteCT_CG"].sum() 
        Meth_value_mean_CHG=float(exons_gene["siteC_CHG"].sum())/exons_gene["siteCT_CHG"].sum()
        Meth_value_mean_CHH=float(exons_gene["siteC_CHH"].sum())/exons_gene["siteCT_CHH"].sum()
        row_exon_mean=[j,Meth_value_mean_CG,Meth_value_mean_CHG,Meth_value_mean_CHH]
        mean_list.append(row_exon_mean)
    AverageMethExon=AverageMethExon.append(mean_list,ignore_index=True)
    AverageMethExon.columns = ["gene_ID","Meth_value_mean_CG","Meth_value_mean_CHG","Meth_value_mean_CHH"]
    print "finish run exon/intron methymean"
    return AverageMethExon




def gtf2bed(gtf):
    localtime = time.asctime( time.localtime(time.time()) )
    print localtime
    #t1 = time.time()
    genePred=gtfToGenePred(gtf)##5.10 seconds
    #t2 = time.time()
    big_bed=AllinfoAnnotation(genePred) ##5 minutes, 48.53 seconds
    #t3 = time.time()
    gb_bed=Create_gb_bed(big_bed) #0.07 seconds
    #t4 = time.time()
    pmt_bed=Create_pmt_bed() #0.32 seconds
    #t5 = time.time()
    exon_bed,intron_bed=CreateExonIntronBed(big_bed)##43 minutes, 9.32 seconds
    localtime = time.asctime( time.localtime(time.time()) )
    print localtime
    #t6 = time.time()
    return big_bed, gb_bed, pmt_bed, exon_bed, intron_bed ##前面都存過了



def Exp(name, exp):
    exp=pd.read_csv(exp,sep="\t",names=["gene_ID","RPKM"])
    exp.to_csv("%s_exp.txt"%(name),sep="\t",index=False,header=False)
    return exp


def WholeProcess(name, cgmap, cutoff, exp, gb_bed, pmt_bed, exon_bed, intron_bed):
    #localtime = time.asctime( time.localtime(time.time()) )
    #print localtime
    t7 = time.time()
    cgmap_cutoff,bwall,bzfile=cgmap2bwbz(name,cgmap,cutoff) ##10 minutes, 11.75 seconds
    #localtime = time.asctime( time.localtime(time.time()) )
    #print localtime
    exp = Exp(name,exp)
    t8 = time.time()
    GeneMeth=RunMethGeneBw(bwall,gb_bed) #46.54 seconds
    GeneMeth.to_csv("%s_GeneMean.txt"%(name),sep="\t",index=False,header=False)
    print "save genemean"
    #localtime = time.asctime( time.localtime(time.time()) )
    #print localtime
    t9 = time.time()
    PromoterMeth=RunMethGeneBw(bwall,pmt_bed) #44.58 seconds
    PromoterMeth.to_csv("%s_PromoterMean.txt"%(name),sep="\t",index=False,header=False)
    print "save PromoterMean"
    #localtime = time.asctime( time.localtime(time.time()) )
    #print localtime
    #exon_bed=pd.read_csv("%sexon_bed"%(name),sep="\t",header=None,dtype={0:str})
    #intron_bed=pd.read_csv("intron_bed"%(name),sep="\t",header=None,dtype={0:str})
    t10 = time.time()
    Exon_meth=RunMethExon(bzfile,exon_bed)##cgmap_cutoff,多一個bzfile>> 5 hour, 17 minutes, 4.17 seconds
    Exon_meth.to_csv("%s_ExonMean.txt"%(name),sep="\t",index=False,header=False) 
    print "save exon mean"
    #localtime = time.asctime( time.localtime(time.time()) )
    #print localtime
    t11 = time.time()
    Intron_meth=RunMethExon(bzfile,intron_bed) ##cgmap_cutoff,t115開始>>> 4 hour, 5 minutes, 51.59 seconds
    Intron_meth.to_csv("%s_IntronMean.txt"%(name),sep="\t",index=False,header=False)
    t12 = time.time()
    print "save intron mean"
    print "save bw %s"%HowManyTime(t7,t8)
    print "save genemean %s"%HowManyTime(t8,t9)
    print "save pmtmean %s"%HowManyTime(t9,t10)
    print "save exonmean %s"%HowManyTime(t10,t11) ##human 2 hour, 45 minutes, 37.04 seconds
    print "save intronmean %s"%HowManyTime(t11,t12)
    #localtime = time.asctime( time.localtime(time.time()) )
    #print localtime


def main():
    parser = get_parser()
    args = parser.parse_args()
    big_bed, gb_bed, pmt_bed, exon_bed, intron_bed = gtf2bed(args.gtffile)
    if args.samplelist == "False":
        WholeProcess(args.samplename,args.cgmapfilename,args.cutoff, args.expressionfile, gb_bed, pmt_bed, exon_bed, intron_bed)
    if args.samplelist != "False":
        Samplelist=pd.read_csv("%s"%(args.samplelist),header=None,sep="\t")
        for i in range(len(Samplelist)):
            WholeProcess(Samplelist.iloc[i,0],Samplelist.iloc[i,1],args.cutoff, Samplelist.iloc[i,2], gb_bed, pmt_bed, exon_bed, intron_bed)




if __name__ == '__main__':
    main()

