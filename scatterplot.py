# coding=UTF-8

'''
import module
'''
import matplotlib
#matplotlib.use('Agg')
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
#sns.set_style("darkgrid")


def get_parser():
    """
    Create a parser and add arguments
    """
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group('Required arguments')
    group1.add_argument("-n","--samplename",type=str,help="the name of the set of data")
    group1.add_argument("-p","--plot",type=str,default="scatter",choices=["scatter","kernel"],help="create scatterplot or kernel density plot, default is 'scatterplot'")
    group1.add_argument("-c","--context",type=str,default="CG",choices=["CG","CHG","CHH"],help="choose the context of methylation, default is 'CG'")
    group1.add_argument("-t","--target",type=str,default="Promoter",choices=["Promoter","Gene_Body","Exon","Intron"],help="choose the target region of methylation, default is 'Promoter'")
    group2 = parser.add_argument_group('Important general arguments')
    group2.add_argument("-cor","--correlation",default="pearson",choices=["False","pearson","spearman"],help="select the type of correlation, default is 'pearson'")
    group2.add_argument("-re0","--skip0",default="False",choices=["True","False"],help="Whether genes with 0 expression value would be included. Default is to include them")
    group2.add_argument("-thrs","--threshold",default="1000",help="Whether skip genes with expression value that is too high, default is to skip genes higher than 1000 expression. If want to include them, please set 'None'")
    group2.add_argument("-xlim","--xlimit",default="False",help="Nemeric zoom in the gene expression value to clearly understand the distribution")
    group2.add_argument("-ylim","--ylimit",default="False",help="Nemeric zoom in the DNA methylation level to clearly understand the distribution")
    group3 = parser.add_argument_group('Graphing arguments')
    #parser.add_argument("-fig","--figuresize",default=(8,6),help="the figure size") ## this can not work
    group3.add_argument("--dotsize",default=30,type=int,help="dotsize, default is 30")
    group3.add_argument("--textsize",default=20,type=int,help="textsize, default is 20")
    group3.add_argument("--ticksize",default=15,type=int,help="ticksize, default is 15")
    group3.add_argument("--labelsize",default=20,type=int,help="labelsize, default is 20")
    group3.add_argument("--titlesize",default=20,type=int,help="titlesize, default is 20")
    ##args=parser.parse_args() ##in main fun
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
		#meth["intronMeth"]=IntronMeth["Meth_value_mean_%s"%(i)]>>cannot because of the number
		methin=pd.merge(meth,IntronMeth[["gene_ID","Meth_value_mean_%s"%(i)]],on=["gene_ID"],how="outer")
		intronname="Meth_value_mean_%s"%(i)
		methin=methin.rename(columns={intronname: 'Intron'})
		methin["Intron"]=methin["Intron"]*100
		res=pd.merge(methin,exp,on=["gene_ID"],how="inner") ##add FPKM add the column name
		All_value["%s"%(i)]=res
	return All_value,name



def scatterplot(All_value,name,plot="scatter",context="CG",target="Gene_Body",corr="pearson",threshold=1000,xlimit=False,ylimit=False,skip0=False,figuresize=(8,6),dotsize=30,textsize=20,ticksize=15,labelsize=20,titlesize=20):
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
				print "select the df", threshold
			if skip0:
				dfcopy=dfcopy.loc[dfcopy["RPKM"]!=0]
				print "remove 0 rpkm"
			dfcopy=dfcopy.dropna(axis=0,how="any") ## dropna=0>>drop row,how="any">>if one na,drop;="all">>only all na,drop
			ar=dfcopy.values ##change df to 2d array
			if plot =="scatter":
				fig=plt.figure(figsize=figuresize)
				ax1=fig.add_subplot(111)
				ax1.scatter(ar[:,0],ar[:,1],s=dotsize)  ##np.log2(ar[:,0]+0.01)
			#if plot == "scatterHis":
			#	ax1=sns.jointplot(x="RPKM",y="%s"%(j),data=dfcopy,color="royalblue")
			if plot == "kernel":
				plt.figure(figsize=figuresize)
				ax1=sns.kdeplot(ar[:,0],ar[:,1],cmap="Blues", shade=True) ##np.log2(ar[:,1]+0.01)
			#if corr=="pearson" and "spearman" and "leastSquares":
			while corr=="pearson" or corr=="spearman" or corr=="leastSquares":  #while corr in ["pearson","spearman","leastSquares"]:
				if corr == "pearson":
					rho, pval = scipy.stats.pearsonr(ar[:,0], ar[:,1]) ##pearson
					print "pearson",rho, pval
					#plt.text(1,1,"R = %2.2f\nP = %2.2f"%(rho,pval),fontsize=10,horizontalalignment='right',verticalalignment='top')
				if corr == "spearman":
					rho, pval = scipy.stats.spearmanr(ar[:,0], ar[:,1],axis=0) ##spearman , axis=0  spearmanr(ar,axis=0)
					print "spearman",rho, pval
					#plt.text(max(ar[:,0]),max(ar[:,1]),"R = %2.2f\nP = %2.2f"%(rho,pval),fontsize=10,horizontalalignment='right',verticalalignment='top')
				if corr == "leastSquares":
					slope, intercept, rho, pval, std_err = scipy.stats.C(ar[:,0], ar[:,1])
					#plt.legend(["R = %.2f\nP = %.2f"%(rho,pval)],loc=1,markerscale=0,fontsize=10)
					#plt.plot(ar[:,0],intercept+slope*ar[:,0],"black")
					## lackkkkkkk a regression eqation
					#plt.text(max(ar[:,0]),96,"R = %2.2f\nP = %2.2f"%(rho,pval),fontsize=10,horizontalalignment='right',verticalalignment='top')
				plt.text(max(ar[:,0]),ylimit*0.95,"R = %-2.3f\nP = %-2.2e"%(rho,pval),fontsize=textsize,horizontalalignment='right',verticalalignment='top')
				break
			ax1.set_title("%s"%(i),fontsize=titlesize)
			highest_value=max(dfcopy["%s"%(j)])###list可能不用改
			if type(ylimit)==int:
				ax1.set_ylim((0, ylimit))
			if ylimit==False:
				ylimit=highest_value*3/2
				if ylimit>100:  ##另一種是此
					ylimit=100
					ax1.set_ylim((0, ylimit))
				else:
					ax1.set_ylim((0, ylimit))
			plt.xticks(fontsize=ticksize)
			plt.yticks(fontsize=ticksize)
			if type(xlimit)==int:
				ax1.set_xlim((0, xlimit))
			ax1.set_ylabel("%s Methylation "%(j)+"(%)",fontsize=labelsize)
			ax1.set_xlabel("Gene Expression Value",fontsize=labelsize) ##li
			plt.savefig('%s_%s_%s_%s.png'%(name,i,j,plot)) ###合在一起我自己存檔會有問題>>但實際上跑不會
			plt.show()




def main():
    parser = get_parser()
    args = parser.parse_args()
    All_value,name=read_meth_exp_file(args.samplename)
    scatterplot(All_value,name,plot=args.plot,context=args.context,target=args.target,corr=args.correlation,threshold=eval(args.threshold),xlimit=eval(args.xlimit),ylimit=eval(args.ylimit),skip0=eval(args.skip0),dotsize=args.dotsize,textsize=args.textsize,ticksize=args.ticksize,labelsize=args.labelsize,titlesize=args.titlesize) #,figuresize=args.figuresize



if __name__ == '__main__':
    main()

