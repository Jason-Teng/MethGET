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
from datetime import datetime
import argparse
import subprocess
import pyBigWig
import time
import os, sys
import scipy
import seaborn as sns
sns.set_style("darkgrid")


def get_parser():
    """
    Create a parser and add arguments
    """
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group('Required arguments')
    group1.add_argument("-n","--samplename",type=str,help="The name of the set of data")
    parser.add_argument("-p","--plot",type=str,default="region",choices=["region","site"],help="choose the format of the metaplot, default is 'region'")
    parser.add_argument("-ma","--metavaluefile",default="False",help="put in the metaplot value file if you have created, default is False")# >>>can choose the gbbed file() directly    
    parser.add_argument("-nb","--numberofgroup",default=5,type=int,help="define how many group to seperate gene expression, default is 5")
    ##parser.add_argument("-c","--context",type=str,default="CG",choices=["CG","CHG","CHH"],help="choose the context of methylation")>>>now is CG, CHG, CHH at the same time，因為畫比較久，所以一次都呈現，並把結果output，其實no的可以一起算，最後畫圖再刪就好
    #parser.add_argument("-t","--target",type=str,default="Promoter",choices=["Gene_Body","Promoter","Exon","Intron"],help="choose the target region of methylation")
    parser.add_argument("-re0","--remove0RPKM",default="False",choices=["True","False"],help="remove genes that their expression value is equal to 0")    
    group3 = parser.add_argument_group('Important general arguments')
    parser.add_argument("-No_bins","--numberofbins",default=60,type=int,help="define the total of bins from upstream to downstream")
    #parser.add_argument("-psn","--posname",default="TSS",choices=["TSS","TES"],type=int,help="choose the specific position of the metapoint plot, default is transcription start site (TSS)")
    #parser.add_argument("-bp","--basepair",default=2000,type=int,help="the basepairs flanking by the chossing position")
    parser.add_argument("-yaxis","--yaxissetting",default="auto",choices=["auto","set100"],help="choose if the ticks of yaxis set on 100 ")    
    group4 = parser.add_argument_group('Graphing arguments')
    #parser.add_argument("-fig","--figuresize",default=(8,6),help="the figure size")
    group4.add_argument("-xtick","--xticksize",default=20,type=int,help="the size of the xticks (upstream, gene, downstream)")
    group4.add_argument("-ytick","--yticksize",default=15,type=int,help="the size of the yticks (methylation level)")
    group4.add_argument("-label","--labelsize",default=20,type=int,help="labelsize")
    group4.add_argument("-title","--titlesize",default=20,type=int,help="titlesize")
    group4.add_argument("-legend","--legendsize",default=15,type=int,help="legendsize")
    ##args=parser.parse_args() ##in main fun
    return parser



def bw_all(name):#(bw_CG,bw_CHG,bw_CHH):
    ## set to dictionary>> so 3 file can use like only one file
    bw_CG=pyBigWig.open("%s_CG.bw"%(name)) ##
    bw_CHG=pyBigWig.open("%s_CHG.bw"%(name))
    bw_CHH=pyBigWig.open("%s_CHH.bw"%(name))
    bwall={}
    bwall["CG"]=bw_CG  
    bwall["CHG"]=bw_CHG
    bwall["CHH"]=bw_CHH
    #bw_CG.isBigWig() ##make sure the type
    #bw_CHG.isBigWig() ##make sure the type
    #bw_CHH.isBigWig() ##make sure the type
    return bwall,name


def readAllFinalValue(name):
    matrixinfo=pd.read_csv("%smetavaluefile.txt"%(name),sep="\t",header=None)
    n=len(matrixinfo[0].unique())-1
    All_final_value=np.zeros(len(matrixinfo[0].unique()))
    All_final_value=list(All_final_value)
    for gr in range(len(matrixinfo[0].unique())):
        All_final_value[gr]={}
        for con in ["CG","CHG","CHH"]:
            info=matrixinfo.loc[(matrixinfo[0]==matrixinfo[0].unique()[gr]) & (matrixinfo[1]==con)]
            info=info.values.tolist()[0][2:]
            All_final_value[gr][con]=info
    return All_final_value



def AverageMethGeneBW(bwall,chr,left,right):
    """
    get target methylation site from CGmap
    """
    Meth_value_mean={}
    for i in ["CG","CHG","CHH"]:
        bw=bwall["%s"%(i)]
        Meth_mean=np.nanmean(bw.values(chr,left-1,right)) ##left -1 can calculate containing the left side
        Meth_value_mean["%s"%(i)]=Meth_mean
    return Meth_value_mean ## 回傳是一個dicttionary



def Metagenebw(bwall,name,chr,sta,end,dir,No_bins):
    """
    calculate average methylation of each bin
    """
    a=No_bins
    if dir=="+":
        length=end-sta #gene length
        upstream=sta-length*0.5#upstream position>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>可以改這個嗎
        downstream=end+length*0.5 #downstream position
        bins=float((downstream - upstream)/a)#bin's length
        All_position=np.zeros(a+1) 
        All_position[0]=upstream 
        for k in range(1,a+1):####只有60個，但因為add第一個0
            All_position[k]=All_position[k-1]+bins
        left=All_position[0]
        right=All_position[a]
    else:
        length=end-sta #gene length
        upstream=end+length*0.5
        downstream=sta-length*0.5
        bins=float((upstream - downstream)/a)#bin's length
        All_position=np.zeros(a+1) 
        All_position[0]=upstream 
        for k in range(1,a+1):####只有60個，但因為add第一個0，總共有61個
            All_position[k]=All_position[k-1]-bins
        left=All_position[a]
        right=All_position[0]
    Meth_value_mean={}
    for con in ["CG","CHG","CHH"]:
        Meth_value_mean["%s"%(con)]=np.zeros(a)
        if dir=="+":
            for i in range(0,a):  ## All pos is not int, its contain Decimal point
                Meth_value_mean["%s"%(con)][i]=AverageMethGeneBW(bwall,chr,int(All_position[i])-1,int(All_position[i+1])-1)["%s"%(con)]  ###All_position[i]-1(because of bw method of access value),All_position[i+1]-1(because did not contain the last value)
        else:
            for i in range(0,a):
                Meth_value_mean["%s"%(con)][i]=AverageMethGeneBW(bwall,chr,int(All_position[i+1]),int(All_position[i]))["%s"%(con)] ##pos[i] value 比 pos[i]+1 還要高
    return Meth_value_mean



def All_value_metageneplot(bwall,name,read_bed,n=5,No_bins=60,remove0RPKM=False):    
    order=["1st","2nd","3rd"]
    for i in range(4,n+1):
        order.append("%dth"%(i))  ##  想要有幾個就有幾個
    groupname=order[:n]### order 沒有加no    
    no_exp=read_bed.loc[read_bed["RPKM"]==0] ##value=0
    other_exp=read_bed.loc[read_bed["RPKM"]!=0] ##vlaue!=0
    sort_other_exp=other_exp.sort_values(by='RPKM')
    sepvalue=[0]  ##區分exp value的那個值，第一個值是0
    groups={} ### save different bed file
    for k in range(n):
        sepvalue.append(sort_other_exp.iloc[(len(sort_other_exp)*(k+1)/n)-1,4]) ###除法是怎麼算的，決定求出位置後減一了
        groups[order[k]]=sort_other_exp.loc[(sort_other_exp["RPKM"]>sepvalue[k]) & (sort_other_exp["RPKM"]<=sepvalue[k+1])]
    """
    create final value
    """
    All_final_value=[]
    if remove0RPKM==False:
        groupname.insert(0,"no")
        no_exp_All_MethLevel={}
        no_exp_Final_value={} ##一次都產生，最後才
        for con in ["CG","CHG","CHH"]:
            no_exp_All_MethLevel["%s"%(con)]=np.zeros((len(no_exp),No_bins))
            for j in range(1,len(no_exp)+1):    
                try:
                    no_exp_All_MethLevel["%s"%(con)][j-1,:]=Metagenebw(bwall,name,no_exp.iloc[j-1,0],no_exp.iloc[j-1,1],no_exp.iloc[j-1,2],no_exp.iloc[j-1,5],No_bins)["%s"%(con)]
                    print no_exp.iloc[j-1,:], j
                except:
                    print no_exp.iloc[j-1,:], "error@no", j
            no_exp_Final_value["%s"%(con)]=np.zeros(No_bins)
            no_exp_All_MethLevel["%s"%(con)]=np.ma.masked_array(no_exp_All_MethLevel["%s"%(con)],np.isnan(no_exp_All_MethLevel["%s"%(con)]))#掩碼數組
            for m in range(0,No_bins):
                no_exp_Final_value["%s"%(con)][m]=no_exp_All_MethLevel["%s"%(con)][:,m].mean()*100
        All_final_value.append(no_exp_Final_value)
    for i in range(n):  
        other_exp_All_MethLevel={}
        other_exp_Final_value={}
        for con in ["CG","CHG","CHH"]:
            other_exp_All_MethLevel["%s"%(con)]=np.zeros((len(groups[order[i]]),No_bins))
            for j in range(1,len(groups[order[i]])+1):  
                try:
                    other_exp_All_MethLevel["%s"%(con)][j-1,:]=Metagenebw(bwall,name,groups[order[i]].iloc[j-1,0],groups[order[i]].iloc[j-1,1],groups[order[i]].iloc[j-1,2],groups[order[i]].iloc[j-1,5],60)["%s"%(con)]
                    print groups[order[i]].iloc[j-1,:], j
                except:
                    print groups[order[i]].iloc[j-1,:], "error@others", j
            other_exp_Final_value["%s"%(con)]=np.zeros(No_bins)
            other_exp_All_MethLevel["%s"%(con)]=np.ma.masked_array(other_exp_All_MethLevel["%s"%(con)],np.isnan(other_exp_All_MethLevel["%s"%(con)]))#掩碼數組
            for m in range(0,No_bins):
                other_exp_Final_value["%s"%(con)][m]=other_exp_All_MethLevel["%s"%(con)][:,m].mean()*100
        All_final_value.append(other_exp_Final_value)
    matrixinfo=pd.DataFrame()
    matrix = []
    for gr in range(len(groupname)):
        for con in ["CG","CHG","CHH"]:
            info = list(All_final_value[gr][con])
            info.insert(0,"%s"%(con))
            info.insert(0,"%s"%(groupname[gr]))
            matrix.append(info)
    matrixinfo = matrixinfo.append(matrix)
    matrixinfo.to_csv("%smetavaluefile.txt"%(name),header=False,index=False,sep="\t")
    return All_final_value
    


def metageneplot(All_final_value,name,n=5,No_bins=60,remove0RPKM=False,yaxis="auto",figuresize=(8,6),xticksize=20,yticksize=15,labelsize=20,titlesize=20,legendsize=15):
    order=["1st","2nd","3rd"]
    for i in range(4,n+1):
        order.append("%dth"%(i))  ##  想要有幾個就有幾個
    groupname=order[:n]### order 沒有加no
    if remove0RPKM==False:
        groupname.insert(0,"no")  ##分群的名稱，要加no expression, groupname是有加no的    
    for con in ["CG","CHG","CHH"]:
        ####max###
        maxall=[]
        for j in range(len(All_final_value)):  ##+1的是no_exp的
            All_final_value[j]
            maxall.append(max(All_final_value[j]["%s"%(con)]))
        highest_value=max(maxall)
        ####numerical_equation####
        #numerical_equation="no (RPKM=0.00)"
        #for j in range(len(order[:n])):
        #   numerical_equation+=", %s (%-3.2f<RPKM<=%-3.2f)"%(order[:n][j],sepvalue[j],sepvalue[j+1])
        ####plotting####
        plt.figure(figsize=figuresize)
        plt.subplot(111)
        ### ylim(set the value of y lim) ###
        ylimit=highest_value*3/2
        if ylimit>100:  ##另一種是此
            plt.ylim((0, 100))
        else:
            plt.ylim((0, ylimit))
        l1=plt.axvline(x=15, color='lightgray',linestyle='--')
        l2=plt.axvline(x=45, color='lightgray',linestyle='--')
        color=cm.rainbow(np.linspace(0.6,0.0,len(All_final_value)))
        color = color.tolist()
        for j in range(len(All_final_value)):
            plt.plot(All_final_value[j]["%s"%(con)],color=color[j],label=groupname[j])
        plt.legend(loc='upper right',fontsize=legendsize)#fontsize="x-small"
        plt.title("%s"%(con),fontsize=titlesize)
        plt.xticks([8, 30, 52],[r'$upstream$', r'$gene$', r'$downstream$'],fontsize=xticksize) ###改
        plt.yticks(fontsize=yticksize)
        plt.ylabel('Methylation Level (%)',fontsize=labelsize)
        plt.savefig('%s_meta_gene_plot_%s_%dexpress.png'%(name,con,n))
        plt.show()






def Metapointbw(bwall,chr,pos,dir,bp,No_bins):
    totalBin=2*No_bins ##bp and No_bins is one side, therefore total bp and bins are double
    if dir=="+":
        upstream=pos-bp#upstream position
        downstream=pos+bp #downstream position
        bins=float(bp/No_bins)#bin's length
        All_position=np.zeros(totalBin+1) 
        All_position[0]=upstream
        for k in range(1,totalBin+1):####只有60個，但因為add第一個0
            All_position[k]=All_position[k-1]+bins
        left=All_position[0] ## >>>不知道需不需要
        right=All_position[totalBin]
    else:
        #length=end-sta #gene length
        upstream=pos+bp
        downstream=pos-bp
        bins=float(bp/No_bins) #bin's length
        All_position=np.zeros(totalBin+1) 
        All_position[0]=downstream
        for k in range(1,totalBin+1): ####只有60個，但因為add第一個0，總共有61個
            All_position[k]=All_position[k-1]-bins
        left=All_position[totalBin]
        right=All_position[0]
    Meth_value_mean={} ## this dict for CG,CHG,CHH, the average meth in that range
    for con in ["CG","CHG","CHH"]: ##loop for three context
        #Meth_site_context=Meth_site_cut[Meth_site_cut["3letter"]==con]
        Meth_value_mean["%s"%(con)]=np.zeros(totalBin)
        if dir=="+":
            for i in range(0,totalBin):  
                #Meth_select_left=Meth_site_context[Meth_site_context.site>=All_position[i]]  ##有包含前一個值
                #Meth_select_right=Meth_select_left[Meth_select_left.site<All_position[i+1]]  ##沒有包含後一個值
                #Meth_value_mean["%s"%(con)][i]=Meth_select_right.loc[:,"Meth_level"].mean()  ##應該是60個值吧，求methylation level那列平均
                Meth_value_mean["%s"%(con)][i]=AverageMethGeneBW(bwall,chr,int(All_position[i])-1,int(All_position[i+1])-1)["%s"%(con)]
        else:
            for i in range(0,totalBin):
                #Meth_select_right=Meth_site_context[Meth_site_context.site<=All_position[i]] #從後面數過來
                #Meth_select_left=Meth_select_right[Meth_select_right.site>All_position[i+1]]
                #Meth_value_mean["%s"%(con)][i]=Meth_select_left.loc[:,"Meth_level"].mean()
                Meth_value_mean["%s"%(con)][i]=AverageMethGeneBW(bwall,chr,int(All_position[i+1]),int(All_position[i]))["%s"%(con)]
    return Meth_value_mean






def All_value_metapointplot(bwall,name,read_bed,n=5,posname="TSS",bp=2000,No_bins=10,remove0RPKM=False):
    totalBin=2*No_bins
    order=["1st","2nd","3rd"]
    for i in range(4,n+1):
        order.append("%dth"%(i))  ##  想要有幾個就有幾個
    groupname=order[:n]### order 沒有加no    
    no_exp=read_bed.loc[read_bed["RPKM"]==0] ##value=0
    other_exp=read_bed.loc[read_bed["RPKM"]!=0] ##vlaue!=0
    sort_other_exp=other_exp.sort_values(by='RPKM')
    sepvalue=[0]  ##區分exp value的那個值，第一個值是0
    groups={} ### save different bed file
    for k in range(n):
        sepvalue.append(sort_other_exp.iloc[(len(sort_other_exp)*(k+1)/n)-1,4]) ###除法是怎麼算的，決定求出位置後減一了
        groups[order[k]]=sort_other_exp.loc[(sort_other_exp["RPKM"]>sepvalue[k]) & (sort_other_exp["RPKM"]<=sepvalue[k+1])]
    """
    create final value
    """
    All_final_value=[]
    if remove0RPKM==False:
        groupname.insert(0,"no")
        no_exp_All_MethLevel={}
        no_exp_Final_value={} ##一次都產生，最後才 
        for con in ["CG","CHG","CHH"]:
            no_exp_All_MethLevel["%s"%(con)]=np.zeros((len(no_exp),totalBin))
            for j in range(1,len(no_exp)+1):    
                if posname=="TSS":
                    pos=no_exp.iloc[j-1,1]
                if posname=="TES":
                    pos=no_exp.iloc[j-1,2]
                #print(pos) #check pos
                try:
                    no_exp_All_MethLevel["%s"%(con)][j-1,:]=Metapointbw(bwall,no_exp.iloc[j-1,0],pos,no_exp.iloc[j-1,5],bp,No_bins)["%s"%(con)] ###problem 應該是哪裡已經滿了
                    print no_exp.iloc[j-1,:], j
                except:
                    print no_exp.iloc[j-1,:], "error@no", j## because of MT
            no_exp_Final_value["%s"%(con)]=np.zeros(totalBin)
            print no_exp_Final_value["%s"%(con)]
            no_exp_All_MethLevel["%s"%(con)]=np.ma.masked_array(no_exp_All_MethLevel["%s"%(con)],np.isnan(no_exp_All_MethLevel["%s"%(con)]))#掩碼數組>>>這個應該不用
            for m in range(0,totalBin):
                no_exp_Final_value["%s"%(con)][m]=no_exp_All_MethLevel["%s"%(con)][:,m].mean()*100#>>這個改nanmean
        All_final_value.append(no_exp_Final_value)
    for i in range(n):  
        other_exp_All_MethLevel={}
        other_exp_Final_value={}
        for con in ["CG","CHG","CHH"]:
            other_exp_All_MethLevel["%s"%(con)]=np.zeros((len(groups[order[i]]),totalBin))
            for j in range(1,len(groups[order[i]])+1):  
                if posname=="TSS":
                    pos=groups[order[i]].iloc[j-1,1]
                if posname=="TES":
                    pos=groups[order[i]].iloc[j-1,2]
                try:
                    other_exp_All_MethLevel["%s"%(con)][j-1,:]=Metapointbw(bwall,groups[order[i]].iloc[j-1,0],pos,groups[order[i]].iloc[j-1,5],bp,No_bins)["%s"%(con)]
                    print other_exp.iloc[j-1,:], j,i
                except:
                    groups[order[i]].iloc[j-1,:], "error@others", j
            other_exp_Final_value["%s"%(con)]=np.zeros(totalBin)
            other_exp_All_MethLevel["%s"%(con)]=np.ma.masked_array(other_exp_All_MethLevel["%s"%(con)],np.isnan(other_exp_All_MethLevel["%s"%(con)]))#掩碼數組
            for m in range(0,totalBin):
                other_exp_Final_value["%s"%(con)][m]=other_exp_All_MethLevel["%s"%(con)][:,m].mean()*100
        All_final_value.append(other_exp_Final_value)
    matrixinfo=pd.DataFrame()
    matrix = []
    for gr in range(len(groupname)):
        for con in ["CG","CHG","CHH"]:
            info = list(All_final_value[gr][con])
            info.insert(0,"%s"%(con))
            info.insert(0,"%s"%(groupname[gr]))
            matrix.append(info)
    matrixinfo = matrixinfo.append(matrix)
    matrixinfo.to_csv("%s_meta%sfile.txt"%(name,posname),header=False,index=False,sep="\t")
    return All_final_value




def metapointplot(All_final_value,posname="TSS",n=5,bp=2000,No_bins=10,remove0RPKM=False,yaxis="auto",figuresize=(8,6),xticksize=20,yticksize=15,labelsize=20,titlesize=20,legendsize=15):
    order=["1st","2nd","3rd"]
    for i in range(4,n+1):
        order.append("%dth"%(i))  ##  想要有幾個就有幾個
    groupname=order[:n]### order 沒有加no
    if remove0RPKM==False:
        groupname.insert(0,"no")  ##分群的名稱，要加no expression, groupname是有加no的  
    for i in ["CG","CHG","CHH"]:
        #### max ####
        maxall=[]
        for j in range(n+1):  ##+1的是no_exp的
            All_final_value[j] ###this will print out
            maxall.append(max(All_final_value[j]["%s"%(i)]))
        highest_value=max(maxall)
        ####numerical_equation####
        #numerical_equation="no (RPKM=0.00)"
        #for j in range(len(order[:n])):
        #   numerical_equation+=", %s (%-3.2f<RPKM<=%-3.2f)"%(order[:n][j],sepvalue[j],sepvalue[j+1])
        ####plotting####
        plt.figure(figsize=figuresize)
        plt.subplot(111)
        ### ylim(set the value of y lim) ###
        ylimit=highest_value*3/2
        if ylimit>100:  ##另一種是此
            plt.ylim((0, 100))
        else:
            plt.ylim((0, ylimit))
        l1=plt.axvline(x=No_bins-0.5, color='lightgray',linestyle='--') ##x=a-0.5>>>from 0 begin>>>why 20個點,畫在9.5的位置，沒錯的
        color=cm.rainbow(np.linspace(0.6,0.0,n+1))
        color = color.tolist()
        for j in range(n+1):
            plt.plot(All_final_value[j]["%s"%(i)],color=color[j],label=groupname[j])
        plt.legend(loc='upper right',fontsize=legendsize)#fontsize="x-small"
        plt.title("%s"%(i),fontsize=titlesize)
        plt.xticks([0,No_bins-0.5,2*No_bins-1],["-%d"%(bp),'%s'%(posname),"+%d"%(bp)],fontsize=xticksize) ###改過的
        plt.yticks(fontsize=yticksize)
        plt.ylabel('Methylation Level (%)',fontsize=labelsize)
        #plt.savefig('%s_meta_%s_plot_%s_%dexpress_%s.png'%(sample,posname,i,n,date))
        plt.show()



def main():
    parser = get_parser()
    args = parser.parse_args()
    exp=pd.read_csv("%s_exp.txt"%(args.samplename),sep="\t",names=["gene_ID","RPKM"])
    read_bed=pd.read_csv('gb.bed', sep='\t',names=["chr","left","right","gene_ID","RPKM","direction"],dtype={"chr":str})
    res=pd.merge(read_bed[["chr","left","right","gene_ID","direction"]],exp,on=["gene_ID"]) ##how='inner'
    res=res[["chr","left","right","gene_ID","RPKM","direction"]]
    bwall,name=bw_all(args.samplename)
    if args.plot == "region":
        if args.metavaluefile=="False":
            All_final_value=All_value_metageneplot(bwall,name,read_bed=res,n=args.numberofgroup,No_bins=args.numberofbins,remove0RPKM=eval(args.remove0RPKM)) #,n=args.numberofgroup,No_bins=args.numberofbins,remove0RPKM=eval(args.remove0RPKM)
        if args.metavaluefile=="True":
            All_final_value=readAllFinalValue(name)
        metageneplot(All_final_value,name,n=args.numberofgroup,No_bins=args.numberofbins,remove0RPKM=eval(args.remove0RPKM),yaxis="auto",figuresize=(8,6),xticksize=args.xticksize,yticksize=args.yticksize,labelsize=args.labelsize,titlesize=args.labelsize,legendsize=args.legendsize)
    if args.plot == "site":
        if args.metavaluefile=="False":
            All_final_value=All_value_metapointplot(bwall,name,read_bed=res,n=args.numberofgroup,No_bins=args.numberofbins,remove0RPKM=eval(args.remove0RPKM)) #,n=args.numberofgroup,No_bins=args.numberofbins,remove0RPKM=eval(args.remove0RPKM)
        if args.metavaluefile=="True":
            All_final_value=readAllFinalValue(name)
        metapointplot(All_final_value,name,n=args.numberofgroup,No_bins=args.numberofbins,remove0RPKM=eval(args.remove0RPKM),yaxis="auto",figuresize=(8,6),xticksize=args.xticksize,yticksize=args.yticksize,labelsize=args.labelsize,titlesize=args.labelsize,legendsize=args.legendsize)




if __name__ == '__main__':
    main()


