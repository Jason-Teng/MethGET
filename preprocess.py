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
import warnings
warnings.filterwarnings('ignore')


def get_parser():
    """
    Create a parser and add arguments
    """
    parser = argparse.ArgumentParser(usage="Usage: preprocess.py {-s <samplelist> | -n <samplename> -f <cgmap> -e <expressionfile>} -g <gtf> [options]") 
    group2 = parser.add_argument_group('For samplelist')
    group2.add_argument("-s","--samplelist",type=str,default="False",help="input sample description file")
    group1 = parser.add_argument_group('For individual data')
    group1.add_argument("-n","--samplename",type=str,help="Setting the name of the set of data. This determines the names of output files for downstream analyses")
    group1.add_argument("-f","--cgmapfilename",type=str,help="input CGmap.gz file")
    group1.add_argument("-e","--expressionfile",type=str,help="input gene expression file")
    group3 = parser.add_argument_group('General options')
    group3.add_argument("-g","--gtffile",type=str,help="input gene annotation file")
    group3.add_argument("-c","--cutoff",type=int,default=5, help="minimum cytosines that are covered by reads")
    group3.add_argument("--TE",default="False",choices=["True","False"], help="The input GTF file is TE GTF. The TE methylation each gene will be calculated in GeneMean.txt. The TE analyses can be conducted if the genomic location is 'Gene_Body'")
    return parser


##########     preprocess methylation    ##########

def cgmap2bwbz(name,cgfile,cutoff):
    """
    cut the depth of cgmap, and save into bw file (CG, CHG, CHH)
    """
    cgmapname = str(cgfile)[:-3]
    subprocess.call("zcat %s>%s"%(cgfile, cgmapname), shell=True)
    cgmap = pd.read_csv(cgmapname, sep='\t',header=None,dtype={0:str}, error_bad_lines=False) 
    cgmap_cutoff=cgmap[cgmap[7]>=cutoff]
    if len(cgmap_cutoff)>0:
    	print "finish depth cutoff"
    else:
    	print "no file"
    cgmap_cutoff.to_csv("%s.CGmap.gz"%(name),header=False,sep="\t",index=False,compression='gzip')
    cgmapc=cgmap_cutoff.copy()
    cgoutfile=name+"_CG.bw"
    chgoutfile=name+"_CHG.bw"
    chhoutfile=name+"_CHH.bw"
    ## the boundary of bw in each chr
    bwheader=[]
    for chr in cgmapc.iloc[:,0].unique():
    	k=cgmapc[cgmapc.iloc[:,0]==chr]
    	maxPos=k.iloc[:,2].max()
    	bwheader.append((chr,int(maxPos)+100000)) 
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
    	CG = k[(k.iloc[:,3] == 'CG') ]
    	CHG = k[(k.iloc[:,3] == 'CHG')]
    	CHH = k[(k.iloc[:,3] == 'CHH')]
        try:
        	bw_CG.addEntries(chromosome,[x-1 for x in CG.ix[:,2].tolist()],values=[float(f) for f in CG.ix[:,5].tolist()], span=1)
        	bw_CHG.addEntries(chromosome,[x-1 for x in CHG.ix[:,2].tolist()],values=[float(f) for f in CHG.ix[:,5].tolist()], span=1)
        	bw_CHH.addEntries(chromosome,[x-1 for x in CHH.ix[:,2].tolist()],values=[float(f) for f in CHH.ix[:,5].tolist()], span=1)
        except:
            pass
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
    return cgmap_cutoff,bwall


##########     preprocess gene(TE) annotation     ##########

def gtfToGenePred(gtf): 
    """
    convert gtf file to genePred file
    """
    genePred="genePred"
    subprocess.call("gtfToGenePred %s genePred -genePredExt"%(gtf),shell=True)
    return genePred



def AllinfoTE(genePred):
    """
    convert TE genepred to TE bed (will put in gb.bed)
    """
    df_a = pd.read_csv("genePred",sep='\t',header=None)
    gb_bed = df_a[[1,3,4,0,2]]
    fill=np.full((len(gb_bed),1),".")
    gb_bed.insert(4,".",fill)
    gb_bed.to_csv("gb.bed",sep="\t",index=False,header=False)
    return gb_bed



def AllinfoAnnotation(genePred):
    """
    convert gene genepred to big bed (having information of gene body, promoter, exon, intron)
    """
    subprocess.call("cat %s|awk \'OFS=\"\t\"{print $12,$1,$2,$3,$4,$5,$5-$4,$8,$9,$10}\'>all_info_annottion"%(genePred),shell=True)
    df_a=pd.read_csv("all_info_annottion",sep='\t',header=None)
    big_bed = pd.DataFrame()
    gene_list=df_a[0].unique()
    for i in range(len(gene_list)):
        gene=df_a[df_a[0]==gene_list[i]] 
        if len(gene)>1:
            ##choose the longest gene
            gene=gene.loc[gene[6]==max(gene[6])]
        if len(gene)>1:
            ##choose the more exons genes
            gene=gene.loc[gene[7]==max(gene[7])]
        onlygene=gene.drop_duplicates(subset=[0], keep='first', inplace=False)
        big_bed=big_bed.append(onlygene) 
    big_bed.to_csv("big_bed",sep="\t",index=False,header=False)
    print("save the information of gene annotation: big_bed")
    print("Total gene number is %d"%(len(big_bed)))
    return big_bed



def Create_gb_bed(big_bed):
    '''
    save gene body bed
    '''
    gb_bed=big_bed[[2,4,5,0,3]]  
    fill=np.full((len(gb_bed),1),".")
    gb_bed.insert(4,".",fill)
    gb_bed.to_csv("gb.bed",sep="\t",index=False,header=False)
    print("save Gene_Body_Bed")
    return gb_bed



def Create_pmt_bed():
    '''
    save PMT bed
    '''
    subprocess.call("cat %s|awk \'OFS=\"\t\" {if  ($6==\"+\") print $1,$2-2000,$2,$4,$5,$6 ;else print $1,$3,$3+2000,$4,$5,$6}\'>pmt.bed"%("gb.bed"),shell=True)
    print("save Promoter_Bed")
    pmt_bed=pd.read_csv("pmt.bed",sep="\t",header=None,dtype={0:str})
    return pmt_bed


 
def CreateExonIntronBed(big_bed):
    """
    save exon and intron bed
    """
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
        # get exons' locations
        for j in range(len(exsta[i])):
            exon_row=[big_bed.iloc[i,2],int(exsta[i][j]),int(exend[i][j]),big_bed.iloc[i,0],big_bed.iloc[i,3]] ###chr,sta,end,gene,dir
            exon_list.append(exon_row)
        # get introns' locations
        for j in range(len(exsta[i])-1):
            intron_row=[big_bed.iloc[i,2],int(exend[i][j]),int(exsta[i][j+1]),big_bed.iloc[i,0],big_bed.iloc[i,3]]
            intron_list.append(intron_row)
    exon_bed = exon_bed.append(exon_list)
    intron_bed = intron_bed.append(intron_list)
    fill_exon=np.full((len(exon_bed),1),".") 
    exon_bed.insert(4,".",fill_exon)
    fill_intron=np.full((len(intron_bed),1),".") 
    intron_bed.insert(4,".",fill_intron)
    exon_bed.to_csv("exon.bed",sep="\t",index=False,header=False)
    print("save exon bed")
    intron_bed.to_csv("intron.bed",sep="\t",index=False,header=False)
    print("save intron bed")
    return exon_bed,intron_bed



def gtf2bed(gtf):
    """
    the main script for gene annotation preprocess
    """
    localtime = time.asctime( time.localtime(time.time()) )
    genePred=gtfToGenePred(gtf)
    big_bed=AllinfoAnnotation(genePred) 
    gb_bed=Create_gb_bed(big_bed) 
    pmt_bed=Create_pmt_bed() 
    exon_bed,intron_bed=CreateExonIntronBed(big_bed)
    localtime = time.asctime( time.localtime(time.time()) )
    return big_bed, gb_bed, pmt_bed, exon_bed, intron_bed 



def TEgtf2bed(gtf):
    """
    the main script for TE annotation preprocess (in gb.bed), create some empty files to avoid error
    """
    genePred=gtfToGenePred(gtf)
    gb_bed = AllinfoTE(genePred)
    big_bed = pd.DataFrame()
    big_bed.to_csv("big.bed")
    pmt_bed = pd.DataFrame()
    pmt_bed.to_csv("pmt.bed")
    exon_bed = pd.DataFrame()
    exon_bed.to_csv("exon.bed")
    intron_bed = pd.DataFrame()
    intron_bed.to_csv("intron.bed") 
    return big_bed, gb_bed, pmt_bed, exon_bed, intron_bed


##########     DNA methylation level in the genomic locations     ###########

def AverageMethGeneBW(bwall,chr,left,right):
    """
    the function to average methylaiton level per location (with CG, cHG, CHH)
    """
    Meth_value_mean={}
    for i in ["CG","CHG","CHH"]:
        bw=bwall["%s"%(i)]
        Meth_mean=np.nanmean(bw.values(chr,left-1,right))
        Meth_value_mean["%s"%(i)]=Meth_mean
    return Meth_value_mean



def RunMethGeneBw(bwall,read_bed):
    """
    calculate all average methylation in each locations from the bed file 
    """
    df1 = pd.DataFrame()
    df1_list=[]
    for j in range(len(read_bed)):
        try:
            Meth_value_mean=AverageMethGeneBW(bwall,str(read_bed.iloc[j,0]),read_bed.iloc[j,1],read_bed.iloc[j,2])
            s1 = list((read_bed.iloc[j,3],Meth_value_mean["CG"],Meth_value_mean["CHG"],Meth_value_mean["CHH"]))
            df1_list.append(s1)
            #print j, s1
        except:
            #print j  ###they may be Pt, ...
            pass
    df1 = df1.append(df1_list, ignore_index=True)
    df1.columns=['gene_ID','Meth_value_mean_CG',"Meth_value_mean_CHG","Meth_value_mean_CHH"]
    return df1


def AverageMethExonBW(bwall,chr,left,right):
    """
    the function to get exon methylation per exon or intron (with CG, cHG, CHH)
    """
    Meth_value_sum={}
    Exon_cytosine_count = {}
    for i in ["CG","CHG","CHH"]:
        bw=bwall["%s"%(i)]
        site = bw.values(chr,left-1,right)
        Meth_sum=np.nansum(site)
        Meth_value_sum["%s"%(i)]=Meth_sum
        Exon_cytosine_count["%s"%(i)] = (np.isnan(site)==False).sum()
    return Meth_value_sum, Exon_cytosine_count


def RunMethExon(bwall,extron_bed):
    """
    acquire exon methylation and intron methylation in each gens from bed
    """
    extron_bed.columns = ["chr","left","right","gene_ID","RPKM","direction"]
    SumCount = pd.DataFrame()
    row_list = []
    # the count and methylation level per exons (introns)
    for i in range(len(extron_bed)):
        try:
            Meth_value_sum,Exon_cytosine_count = AverageMethExonBW(bwall,extron_bed.iloc[i,0],extron_bed.iloc[i,1],extron_bed.iloc[i,2])
            row_value = [extron_bed.iloc[i,3],Meth_value_sum["CG"],Exon_cytosine_count["CG"],Meth_value_sum["CHG"],Exon_cytosine_count["CHG"],Meth_value_sum["CHH"],Exon_cytosine_count["CHH"]]
            row_list.append(row_value)
        except RuntimeError:
            pass
    SumCount = SumCount.append(row_list)
    SumCount.columns = ['gene_ID','sum_CG',"count_CG","sum_CHG","count_CHG","sum_CHH","count_CHH"]
    AverageMethExon = pd.DataFrame()
    mean_list = [] #for final value
    # the exon (intron) methylation per gene
    for j in SumCount["gene_ID"].unique():
        exons_gene=SumCount.loc[SumCount["gene_ID"]==j] ## choose all exons for that gene
        Meth_value_mean_CG = np.array(float(exons_gene["sum_CG"].sum()))/np.array(exons_gene["count_CG"].sum()) 
        Meth_value_mean_CHG = np.array(float(exons_gene["sum_CHG"].sum()))/np.array(exons_gene["count_CHG"].sum())
        Meth_value_mean_CHH = np.array(float(exons_gene["sum_CHH"].sum()))/np.array(exons_gene["count_CHH"].sum())
        row_exon_mean=[j,Meth_value_mean_CG,Meth_value_mean_CHG,Meth_value_mean_CHH]
        mean_list.append(row_exon_mean)
    AverageMethExon=AverageMethExon.append(mean_list,ignore_index=True)
    AverageMethExon.columns = ["gene_ID","Meth_value_mean_CG","Meth_value_mean_CHG","Meth_value_mean_CHH"]
    return AverageMethExon



def Exp(name, exp):
    """
    preprocess expression (just change the name of file)
    """
    exp=pd.read_csv(exp,sep="\t",names=["gene_ID","RPKM"])
    exp.to_csv("%s_exp.txt"%(name),sep="\t",index=False,header=False)
    return exp



def WholeProcess(name, cgmap, cutoff, exp, gb_bed, pmt_bed, exon_bed, intron_bed):
    """
    the main script for DNA methylation preprocess (contexts (CG, CHG, CHH) and genomic locations (gene body, promoter, exon, and intron))
    """
    t7 = time.time()
    cgmap_cutoff,bwall=cgmap2bwbz(name,cgmap,cutoff) 
    exp = Exp(name,exp)
    t8 = time.time()
    print "save bw %s"%HowManyTime(t7,t8)
    GeneMeth=RunMethGeneBw(bwall,gb_bed)
    GeneMeth.to_csv("%s_GeneMean.txt"%(name),sep="\t",index=False,header=False)
    t9 = time.time()
    print "save genemean %s"%HowManyTime(t8,t9)
    try:
        PromoterMeth=RunMethGeneBw(bwall,pmt_bed)
        PromoterMeth.to_csv("%s_PromoterMean.txt"%(name),sep="\t",index=False,header=False)
        t10 = time.time()
        print "save pmtmean %s"%HowManyTime(t9,t10)
        Exon_meth=RunMethExon(bwall, exon_bed)
        Exon_meth.to_csv("%s_ExonMean.txt"%(name),sep="\t",index=False,header=False) 
        t11 = time.time()
        print "save exonmean %s"%HowManyTime(t10,t11) 
        Intron_meth=RunMethExon(bwall, intron_bed) 
        Intron_meth.to_csv("%s_IntronMean.txt"%(name),sep="\t",index=False,header=False)
        t12 = time.time()
        print "save intronmean %s"%HowManyTime(t11,t12)
    except (ValueError), e:
        print(e)
        print(traceback.format_exc())
        print("TE gtf only process the TE methylation in 'gene body'")



def HowManyTime(tbegin,tend):
    """
    to calculate the time to evaluate the speed
    """
    tTotal=tend-tbegin
    tsec=tTotal%60
    ttolmin=tTotal//60
    thour=ttolmin//60
    tmin=ttolmin%60
    suretime="running time is %d hour, %d minutes, %.2f seconds"%(thour,tmin,tsec)
    return suretime



def main():
    script_start = time.time() # start time
    parser = get_parser()
    args = parser.parse_args()
    if args.TE == "False":
        big_bed, gb_bed, pmt_bed, exon_bed, intron_bed = gtf2bed(args.gtffile) # gene GTF preprocess
    if args.TE == "True":
        big_bed, gb_bed, pmt_bed, exon_bed, intron_bed = TEgtf2bed(args.gtffile) # TE GTF preprocess
    if args.samplelist == "False":
        WholeProcess(args.samplename,args.cgmapfilename,args.cutoff, args.expressionfile, gb_bed, pmt_bed, exon_bed, intron_bed) # single methylome preprocess
    ## preprocess from the info in samplelist
    if args.samplelist != "False":
        Samplelist=pd.read_csv("%s"%(args.samplelist),header=None,sep="\t") # read the samplelist file
        for i in range(len(Samplelist)):
            WholeProcess(Samplelist.iloc[i,0],Samplelist.iloc[i,1],args.cutoff, Samplelist.iloc[i,2], gb_bed, pmt_bed, exon_bed, intron_bed) # preprocess each methylome
    script_end = time.time() # end time
    print "preprocess time %s"%HowManyTime(script_start,script_end) # total time for preprocessing




if __name__ == '__main__':
    main()

