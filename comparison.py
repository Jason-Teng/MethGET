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
import warnings
sns.set_style("darkgrid")

warnings.filterwarnings('ignore')

def get_parser():
    """
    Create a parser and add arguments
    """
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group('Required arguments')
    group1.add_argument("-s","--samplelist",default = "samplelist.txt",type=str,help="put in the sample description file")
    group1.add_argument("-c","--context",type=str,default="CG",choices=["CG","CHG","CHH"],help="choose the context of methylation, default is CG")
    group1.add_argument("-t","--target",type=str,default="Promoter",choices=["Promoter","Gene_Body","Exon","Intron","all"],help="choose the target region of methylation, default is Promoter")
    group1.add_argument("-mthr","--meththreshold",default="auto",help="set cutoff of differential methylated genes. default 'auto' uses Δ methylation, CG:10, CHG:1, CHH:1")
    group1.add_argument("-ethr","--expthreshold",default=1.0,type=float,help="set cutoff to identify genes that have expression changes, default uses Δ expression log2FC is 1")# 現在是倍數了
    group2 = parser.add_argument_group('Important general arguments')
    group2.add_argument("-p","--plot",type=str,default="scatter",choices=["scatter","kernel"],help="create scatterplot or kernel density plot, default is scatter")
    group2.add_argument("-cor","--correlation",default="pearson",choices=["False","pearson","spearman"],help="select the type of correlation, default is pearson")
    group2.add_argument("--shownumber",default="False",choices=["False","True"],help="whether to show the number of significant genes, default False is not show")
    group2.add_argument("--methmin",default=None,help="minimum Δ methylation for x-axis")
    group2.add_argument("--methmax",default=None,help="maximum Δ methylation for x-axis")
    group2.add_argument("--expmin",default=None,help="minimum Δ gene expression for y-axis")
    group2.add_argument("--expmax",default=None,help="maximum Δ gene expression for y-axis")
    group3 = parser.add_argument_group('Graphing arguments')
    group3.add_argument("--dotsize",default=20,type=int,help="dotsize")
    group3.add_argument("--textsize",default=25,type=int,help="textsize")
    group3.add_argument("--ticksize",default=15,type=int,help="ticksize")
    group3.add_argument("--labelsize",default=25,type=int,help="labelsize")
    group3.add_argument("--titlesize",default=25,type=int,help="titlesize")
    group3.add_argument("--fontsize",default=1.2,type=int,help="fontsize")
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
    return All_value



def compare_value(Samplelist="samplelist.txt", context="CG", target="Gene_Body", add=1, Methcutoff="auto", Expcutoff=2):
    Samplelist=pd.read_csv(Samplelist,header=None,sep="\t") 
    con=[]
    con.append(context)
    tar=[]
    tar.append(target)
    value={}
    for i in con: 
        for j in tar:
            for k in Samplelist[0]: 
                value[k]=read_meth_exp_file(k)  
            meanvalue=value[k][i][["gene_ID"]]  
            everyvalue=value[k][i][["gene_ID"]]
            treatment=Samplelist[3].unique()
            ##create meanvalue for comparison
            for trt in treatment: 
                formethaverage=pd.DataFrame()
                forexpaverage=pd.DataFrame() 
                Samples=Samplelist.loc[Samplelist[3]==trt]
                for k in range(len(Samples)):
                    formethaverage[Samples.iloc[k,0]]=value[Samples.iloc[k,0]][i][j] 
                    everyvalue["%smeth"%(Samples.iloc[k,0])]=value[Samples.iloc[k,0]][i][j]
                    forexpaverage[Samples.iloc[k,0]]=value[Samples.iloc[k,0]][i]["RPKM"]
                    everyvalue["%sexp"%(Samples.iloc[k,0])]=value[Samples.iloc[k,0]][i]["RPKM"]
                meanvalue["%smethmean"%(trt)]=formethaverage.mean(1).copy()
                meanvalue["%sexpmean"%(trt)]=forexpaverage.mean(1).copy()+add
            dropres=meanvalue.dropna(axis=0,how="any") 
            dropres["Δmethylation"]= dropres["%smethmean"%(treatment[0])].copy()-dropres["%smethmean"%(treatment[1])].copy()
            dropres["Δgene expression"]=dropres["%sexpmean"%(treatment[0])].copy()/dropres["%sexpmean"%(treatment[1])].copy()
            dropres["Δgene expression"]=dropres["Δgene expression"].apply(np.log2)
            deltaMethmean = dropres["Δmethylation"].mean()
            deltaExpmean = dropres["Δgene expression"].mean()
            dropres["var"] = (dropres["Δmethylation"]-deltaMethmean)**2 + (dropres["Δgene expression"]-deltaExpmean)**2
            variation = dropres["var"].sum()/(len(dropres)-2)
            sqrtvar = np.sqrt(variation)
            if Methcutoff == "auto":
                if i == "CG":
                    Methcutoff = float(10)
                if i =="CHG":
                    Methcutoff = float(1)
                if i == "CHH":
                    Methcutoff = float(1) 
            else:
                Methcutoff = float(Methcutoff)
            #log2Expcutoff=np.log2(Expcutoff)
            log2Expcutoff = Expcutoff
            selectMethExp_1 = dropres.loc[(dropres["Δmethylation"]>Methcutoff) & (dropres["Δgene expression"]>log2Expcutoff)]
            selectMethExp_2 = dropres.loc[(dropres["Δmethylation"]<Methcutoff*(-1)) & (dropres["Δgene expression"]>log2Expcutoff)]
            selectMethExp_3 = dropres.loc[(dropres["Δmethylation"]<Methcutoff*(-1)) & (dropres["Δgene expression"]<log2Expcutoff*(-1))]
            selectMethExp_4 = dropres.loc[(dropres["Δmethylation"]>Methcutoff) & (dropres["Δgene expression"]<log2Expcutoff*(-1))]
            selectMeth=dropres.loc[(abs(dropres["Δmethylation"])>Methcutoff)]
            selectMeth.reset_index(inplace=True) 
            selectMethExp=selectMeth.loc[(abs(selectMeth["Δgene expression"])>log2Expcutoff)] 
    return Samplelist,dropres,selectMeth,selectMethExp, selectMethExp_1, selectMethExp_2, selectMethExp_3, selectMethExp_4, deltaMethmean, deltaExpmean, variation



def CompareScatterplot(dropres,selectMethExp,selectMethExp_1,selectMethExp_2,selectMethExp_3,selectMethExp_4, variation, plot="scatter",context="CG",target="Gene_Body",corr="pearson",shownumber=True,methmin=None,methmax=None,expmin=None,expmax=None,figuresize=(10,8),dotsize=15,textsize=20,ticksize=15,labelsize=20,titlesize=20):
    con=[]
    con.append(context)
    tar=[]
    tar.append(target)
    value={}
    for i in con: ##["CG","CHG","CHH"]
        for j in tar: ##["Gene_Body","Promoter","Exon","Intron"]    
            if plot == "scatter":
                fig=plt.figure(figsize=figuresize)
                ax1=fig.add_subplot(111)
                ax1.scatter(list(dropres["Δmethylation"]),list(dropres["Δgene expression"]),s=dotsize)
                ##add select plots
                ax1.scatter(list(selectMethExp["Δmethylation"]),list(selectMethExp["Δgene expression"]),s=dotsize,c="red")             
            a=max(list(dropres["Δmethylation"]))
            b=max(list(dropres["Δgene expression"]))
            c=min(list(dropres["Δmethylation"]))
            d=min(list(dropres["Δgene expression"]))
            if methmax!=None:
                methmax = float(methmax)
                a = methmax*0.95
            if expmax!=None:
                expmax = float(expmax)
                b = expmax*0.95
            if methmin!=None:
                methmin = float(methmin)
                c = methmin*0.95
            if expmin!=None:
                expmin = float(expmin)
                d = expmin*0.95
            if plot == "kernel":
                plt.figure(figsize=figuresize)
                ax1=sns.kdeplot(list(dropres["Δmethylation"]),list(dropres["Δgene expression"]),cmap="Blues", shade=True)
            while corr=="pearson" or corr=="spearman" or corr=="leastSquares":
                if corr == "pearson":
                    rho, pval = scipy.stats.pearsonr(dropres["Δmethylation"], dropres["Δgene expression"])
                if corr == "spearman":
                    rho, pval = scipy.stats.spearmanr(dropres["Δmethylation"], dropres["Δgene expression"])
                if corr == "leastSquares":
                    slope, intercept, rho, pval, std_err = scipy.stats.C(ar[:,0], ar[:,1])
                plt.text(a,b,"R = %-2.3f\nP = %-2.2e"%(rho,pval),fontsize=textsize,horizontalalignment='right',verticalalignment='top')
                break
            ax1.set_title("%s %s"%(i,j),fontsize=titlesize)
            if shownumber == "True":
                plt.text(a,b,"%d"%(len(selectMethExp_1)),fontsize=textsize,horizontalalignment='right',verticalalignment='top')
                plt.text(c,b,"%d"%(len(selectMethExp_2)),fontsize=textsize,horizontalalignment='left',verticalalignment='top')
                plt.text(c,d,"%d"%(len(selectMethExp_3)),fontsize=textsize,horizontalalignment='left',verticalalignment='bottom')
                plt.text(a,d,"%d"%(len(selectMethExp_4)),fontsize=textsize,horizontalalignment='right',verticalalignment='bottom')
            plt.xlabel(r'$\Delta\ methylation$'+" (%)",fontsize=labelsize,)
            plt.ylabel(r'$\Delta\ gene\ expression\ (log2\ fold\ change)$',fontsize=labelsize)
            plt.xticks(fontsize=ticksize)
            plt.yticks(fontsize=ticksize)
            plt.xlim(methmin,methmax)
            plt.ylim(expmin,expmax)
            plt.axhline(0,color="gray",linestyle="--") 
            plt.axvline(0,color="gray",linestyle="--") 
            deltascatterplotname = 'comparison%s_plot.png'%(plot)
            plt.savefig(i + "_" + j +"_"+ deltascatterplotname, dpi = 300)
            deltascattertablename = "chosen_genes.txt" 
            try:
                ratio = float(len(selectMethExp_1)*len(selectMethExp_3))/(len(selectMethExp_2)*len(selectMethExp_4))
                deltainfo = "#quardrant_1: %d, #quardrant_2: %d, #quardrant_3: %d, #quardrant_4: %d, ovaerall variation: %f, ratio (+/-): %f."%(len(selectMethExp_1), len(selectMethExp_2), len(selectMethExp_3), len(selectMethExp_4), variation,ratio)
            except:
                ratio = "division by zero"
                deltainfo = "#quardrant_1: %d, #quardrant_2: %d, #quardrant_3: %d, #quardrant_4: %d, ovaerall variation: %f, ratio (+/-): %s."%(len(selectMethExp_1), len(selectMethExp_2), len(selectMethExp_3), len(selectMethExp_4), variation,ratio)
            print("information for %s, %s"%(context, target), deltainfo)
            df = selectMethExp.iloc[:,1:-1]        
            df.to_csv(i + "_" + j + "_" + deltascattertablename, sep="\t",index=False, header=True) #sep="\t"
            #plt.show()
            plt.close()
    return deltascatterplotname, deltascattertablename



def main():
    parser = get_parser()
    args = parser.parse_args()
    meth_con = [args.context]
    meth_tar = [args.target]
    if args.context == "all":
        meth_con = ["CG", "CHG", "CHH"]
    if args.target == "all":
        meth_tar = ["Promoter","Gene_Body","Exon","Intron"]
    for con in meth_con:
        for tar in meth_tar:
            corr = args.correlation
            if args.shownumber==True:
                corr = False
            Samplelist,dropres,selectMeth,selectMethExp, selectMethExp_1, selectMethExp_2, selectMethExp_3, selectMethExp_4, deltaMethmean, deltaExpmean, variation = compare_value(Samplelist=args.samplelist,context=con, target=tar,add=1.0, Methcutoff=args.meththreshold,Expcutoff= args.expthreshold)
            deltascatterplotname, deltascattertablename = CompareScatterplot(dropres, selectMethExp, selectMethExp_1, selectMethExp_2, selectMethExp_3, selectMethExp_4, variation, plot=args.plot, context=con,target=tar, corr=corr, shownumber=args.shownumber,methmin=args.methmin,methmax=args.methmax,expmin=args.expmin,expmax=args.expmax,figuresize=(10,8),dotsize=args.dotsize,textsize=args.textsize,ticksize=args.ticksize,labelsize=args.labelsize,titlesize=args.titlesize)





if __name__ == '__main__':
    main()

