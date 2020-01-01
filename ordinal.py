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
    group1.add_argument("-n","--samplename", type=str, help="the name of the set of data")
    group1.add_argument("-c","--context",type=str,default="all",choices=["CG","CHG","CHH","all"],help="choose the context of methylation, default 'all' is to choose them all")
    group1.add_argument("-t","--target",type=str,default="all",choices=["Gene_Body","Promoter","Exon","Intron","all"],help="choose the target region of methylation, default 'all' is to choose them all")
    group3 = parser.add_argument_group('Important general arguments')  
    group3.add_argument("-re0","--skip0",default="False",choices=["True","False"],help="whether genes with 0 expression value would be included. Default 'False' is to include them")
    group3.add_argument("-thrs","--threshold",default="2000",help="whether to skip genes with expression value that is too high, default is to skip genes higher than 2000. If want to include them, please set 'None'")
    group2 = parser.add_argument_group('Chart visulaization arguments')
    group2.add_argument("-shsca","--showscatterplot",default="True",choices=["True","False"],help="whether to show the scatterplot, default is to show")
    group2.add_argument("-line","--smoothline",default="True",choices=["True","False"],help="whether to show the fitting curves, default is to show")
    group2.add_argument("-smoo","--smooth_N",default=500,type=int,help="set the number of ticks to average when drawing the fitting curve, default is 500")
    group2.add_argument("-ylim","--ylimit",default="False",help="numeric zoom in the DNA methylation level to clearly understand the distribution")
    group4 = parser.add_argument_group('Graphing arguments')
    group4.add_argument("--dotsize",default=20,type=int,help="dotsize, default is 20")
    group4.add_argument("--textsize",default=25,type=int,help="textsize, default is 25")
    group4.add_argument("--ticksize",default=15,type=int,help="ticksize, default is 15")
    group4.add_argument("--labelsize",default=25,type=int,help="labelsize,, default is 25")
    group4.add_argument("--titlesize",default=25,type=int,help="titlesize, default is 25")
    group4.add_argument("--legendsize",default=20,type=int,help="legendsize, default is 20")
    return parser


def read_meth_exp_file(name):
    exp=pd.read_csv("%s_exp.txt"%(name),sep="\t",names=["gene_ID","RPKM"])
    GeneMeth=pd.read_csv("%s_GeneMean.txt"%(name),sep="\t",names=["gene_ID","Meth_value_mean_CG","Meth_value_mean_CHG","Meth_value_mean_CHH"])
    PromoterMeth=pd.read_csv("%s_PromoterMean.txt"%(name),sep="\t",names=["gene_ID","Meth_value_mean_CG","Meth_value_mean_CHG","Meth_value_mean_CHH"])
    ExonMeth=pd.read_csv("%s_ExonMean.txt"%(name),sep="\t",names=["gene_ID","Meth_value_mean_CG","Meth_value_mean_CHG","Meth_value_mean_CHH"])
    IntronMeth=pd.read_csv("%s_IntronMean.txt"%(name),sep="\t",names=["gene_ID","Meth_value_mean_CG","Meth_value_mean_CHG","Meth_value_mean_CHH"])
    All_value={}
    for i in ["CG","CHG","CHH"]:
        meth=pd.DataFrame()
        meth["gene_ID"]=GeneMeth["gene_ID"]
        meth["Gene_Body"]=GeneMeth["Meth_value_mean_%s"%(i)]*100
        meth["Promoter"]=PromoterMeth["Meth_value_mean_%s"%(i)]*100
        meth["Exon"]=ExonMeth["Meth_value_mean_%s"%(i)]*100
        methin=pd.merge(meth,IntronMeth[["gene_ID","Meth_value_mean_%s"%(i)]],on=["gene_ID"],how="outer")
        intronname="Meth_value_mean_%s"%(i)
        methin=methin.rename(columns={intronname: 'Intron'})
        methin["Intron"]=methin["Intron"]*100
        res=pd.merge(methin,exp,on=["gene_ID"],how="inner") 
        All_value["%s"%(i)]=res
    return All_value,name


def szooth(x,bin):
    afsmoo=[]
    for i in range(len(x)-bin+1):
        n=0
        for j in range(bin):
            n+=x[i+j]
        afsmoo.append(float(n)/float(bin))
    return afsmoo



def rankscatter(All_value,name,plot="scatter",context="CG",target="Gene Body",threshold=2000,skip0=False,RPKM0="on",showscatterplot=True,ylimit=False,smoothline=True,smooth_N=500,figuresize=(11,8),dotsize=30,textsize=20,ticksize=15,labelsize=25,titlesize=25,legendsize=20):
    con=[]
    con.append(context)
    tar=[]
    tar.append(target)
    for i in con:##["CG","CHG","CHH"]
        for j in tar:##["Gene Body","Promoter","Exon","Intron"]
            sort_all=All_value["%s"%(i)].sort_values(by="RPKM")
            #y=sort_all["pmtMeth"]
            pdcopy=sort_all[["%s"%(j),"RPKM"]].copy()
            if type(threshold)==int:
                pdcopy=pdcopy.loc[pdcopy["RPKM"]<threshold]
            if skip0:
                pdcopy=pdcopy.loc[pdcopy["RPKM"]!=0]
            pdcopy=pdcopy.dropna(axis=0,how="any") 
            y=list(pdcopy["%s"%(j)])
            x=np.arange(len(y))
            if plot=="scatter":
                fig=plt.figure(figsize=figuresize)
                ax1=fig.add_subplot(111)
                plt.setp(ax1.get_yticklabels(), fontsize=ticksize)
                ax1.scatter(x,y,color="lightskyblue",s=1,marker="o") 
                if smoothline:
                    yhat = szooth(y,smooth_N)
                    for l in range(smooth_N/2):
                        yhat.insert(0,np.nan)
                    ln1=ax1.plot(yhat,color="blue",label="methylation")
                if smoothline:
                    ax2=ax1.twinx()
                    ln2=ax2.plot(list(pdcopy["RPKM"]),c="r",label="expression")
                    ax2.set_ylabel("Gene Expression Value",size=labelsize)
                    plt.setp(ax2.get_yticklabels(), fontsize=ticksize)
                    lns=ln1+ln2
                    labs = [l.get_label() for l in lns]
                    ax1.legend(lns, labs, loc=1,bbox_to_anchor=(0.92,0.95),fontsize=legendsize)
            if skip0 == False:
                minexp=list(pdcopy.loc[pdcopy["RPKM"]!=0]["RPKM"])[0]
                pos=list(pdcopy["RPKM"]).index(minexp)
                plt.axvline(pos,color="gray",linestyle="--")
            ax1.set_ylabel("%s Methylation "%(j)+"(%)",fontsize=labelsize)
            ax1.set_xlabel("Ranked Genes by Expression Level",size=labelsize)
            try:
                if type(ylimit)==int:
                    ax1.set_ylim((0, ylimit))  
                del ylimit
            except:
                pass        
            plt.xticks([],fontsize=ticksize)
            ax1.set_title("%s"%(i),fontsize=titlesize)
            plt.savefig('%s_%s_%s_smooth_graph_diagram.png'%(name,i,j), dpi=300)
            print("complete %s on %s"%(context,target))
            #plt.show()



def main():
    parser = get_parser()
    args = parser.parse_args()
    All_value,name=read_meth_exp_file(args.samplename)
    meth_con = [args.context]
    meth_tar = [args.target]
    if args.context == "all":
        meth_con = ["CG", "CHG", "CHH"]
    if args.target == "all":
        meth_tar = ["Promoter","Gene_Body","Exon","Intron"]
    for con in meth_con:
        for tar in meth_tar:
            rankscatter(All_value,name,plot="scatter",context=con,target=tar,threshold=eval(args.threshold),skip0=eval(args.skip0),RPKM0="on",showscatterplot=eval(args.showscatterplot),ylimit=eval(args.ylimit),smoothline=eval(args.smoothline),smooth_N=args.smooth_N,figuresize=(11,8),dotsize=args.dotsize,textsize=args.textsize,ticksize=args.ticksize,labelsize=args.labelsize,titlesize=args.titlesize,legendsize=args.legendsize)
    


if __name__ == '__main__':
    main()
