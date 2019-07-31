# coding=UTF-8

'''
import module
'''
import matplotlib
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
    group1.add_argument("-n","--samplename",type=str,help="the name of the set of data")
    group1.add_argument("-p","--plot",type=str,default="scatter",choices=["scatter","kernel"],help="create scatterplot or kernel density plot, default is scatter")
    group1.add_argument("-c","--context",type=str,default="all",choices=["CG","CHG","CHH", "all"],help="choose the context of methylation, default 'all' is to choose them all")
    group1.add_argument("-t","--target",type=str,default="all",choices=["Promoter","Gene_Body","Exon","Intron","all"],help="choose the target region of methylation, default 'all' is to choose them all")
    group2 = parser.add_argument_group('Important general arguments')
    group2.add_argument("-cor","--correlation",default="pearson",choices=["False","pearson","spearman"],help="select the type of correlation, default is pearson")
    group2.add_argument("-re0","--skip0",default="False",choices=["True","False"],help="whether genes with 0 expression value would be included. Default 'False' is to include them")
    group2.add_argument("-thrs","--threshold",default="2000",help="whether skip genes with expression value that is too high, default is to skip genes higher than 2000. If want to include them, please set 'None'")
    group2.add_argument("-xlim","--xlimit",default="False",help="numeric zoom in the gene expression value to clearly understand the distribution")
    group2.add_argument("-ylim","--ylimit",default="False",help="numeric zoom in the DNA methylation level to clearly understand the distribution")
    group3 = parser.add_argument_group('Graphing arguments')
    group3.add_argument("--dotsize",default=30,type=int,help="dotsize, default is 30")
    group3.add_argument("--textsize",default=20,type=int,help="textsize, default is 20")
    group3.add_argument("--ticksize",default=15,type=int,help="ticksize, default is 15")
    group3.add_argument("--labelsize",default=20,type=int,help="labelsize, default is 20")
    group3.add_argument("--titlesize",default=20,type=int,help="titlesize, default is 20")
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



def scatterplot(All_value,name,plot="scatter",context="CG",target="Gene_Body",corr="pearson",threshold=2000,xlimit=False,ylimit=False,skip0=False,figuresize=(8,6),dotsize=30,textsize=20,ticksize=15,labelsize=20,titlesize=20):
	con=[]
	con.append(context)
	tar=[]
	tar.append(target)
	for i in con: ##["CG","CHG","CHH"]
		for j in tar: ##["Gene_Body","Promoter","Exon","Intron"]
			data=All_value["%s"%(i)]
			dfcopy=data[["RPKM","%s"%(j)]].copy()
			if type(threshold)==int:
				dfcopy=dfcopy.loc[dfcopy["RPKM"]<threshold]
			if skip0:
				dfcopy=dfcopy.loc[dfcopy["RPKM"]!=0]
			dfcopy=dfcopy.dropna(axis=0,how="any") 
			ar=dfcopy.values 
			if plot =="scatter":
				fig=plt.figure(figsize=figuresize)
				ax1=fig.add_subplot(111)
				ax1.scatter(ar[:,0],ar[:,1],s=dotsize)  
			if plot == "kernel":
				plt.figure(figsize=figuresize)
				ax1=sns.kdeplot(ar[:,0],ar[:,1],cmap="Blues", shade=True) 
			ax1.set_title("%s"%(i),fontsize=titlesize)
			highest_value=max(dfcopy["%s"%(j)])
			if type(ylimit)==int:
				ax1.set_ylim((0, ylimit))
			if ylimit==False:
				ylimit=highest_value*3/2
				if ylimit > 100:  
					ylimit = 100
					ax1.set_ylim((0, ylimit))
				else:
					ax1.set_ylim((0, ylimit))			
			while corr=="pearson" or corr=="spearman" or corr=="leastSquares": 
				if corr == "pearson":
					rho, pval = scipy.stats.pearsonr(ar[:,0], ar[:,1])
				if corr == "spearman":
					rho, pval = scipy.stats.spearmanr(ar[:,0], ar[:,1],axis=0) 
				plt.text(max(ar[:,0]),ylimit*0.95,"R = %-2.3f\nP = %-2.2e"%(rho,pval),fontsize=textsize,horizontalalignment='right',verticalalignment='top')
				print("complete %s on %s, correlation is: R = %-2.3f, P = %-2.2e"%(context,target,rho,pval))
				break
			plt.xticks(fontsize=ticksize)
			plt.yticks(fontsize=ticksize)
			if type(xlimit)==int:
				ax1.set_xlim((0, xlimit))
			ax1.set_ylabel("%s Methylation "%(j)+"(%)",fontsize=labelsize)
			ax1.set_xlabel("Gene Expression Value",fontsize=labelsize) 
			plt.savefig('%s_%s_%s_%s.png'%(name,i,j,plot), dpi=300) 
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
            scatterplot(All_value,name,plot=args.plot,context=con,target=tar,corr=args.correlation,threshold=eval(args.threshold),xlimit=eval(args.xlimit),ylimit=eval(args.ylimit),skip0=eval(args.skip0),dotsize=args.dotsize,textsize=args.textsize,ticksize=args.ticksize,labelsize=args.labelsize,titlesize=args.titlesize) 



if __name__ == '__main__':
    main()

