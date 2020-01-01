# coding=UTF-8

'''
import module
'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
import shlex, subprocess
from subprocess import Popen, PIPE
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
	group1.add_argument("-n","--samplename",type=str,help="the name of the set of data")
	group1.add_argument("-p","--plot",type=str,default="boxplot",choices=["boxplot","violinplot"],help="create boxplot or violinplot, default is boxplot")
	group1.add_argument("-c","--context",type=str,default="all",choices=["CG","CHG","CHH","all"],help="choose the context of methylation, default 'all' is to choose them all")
	group1.add_argument("-t","--target",type=str,default="all",choices=["Gene_Body","Promoter","Exon","Intron","all"],help="choose the target region of methylation, default 'all' is to choose them all")
	group1.add_argument("-nb","--numberofgroup",default=5,type=int,help="define how many group to seperate gene expression, default is 5")
	group3 = parser.add_argument_group('Important general arguments')
	group3.add_argument("-re0","--skip0",default="False",choices=["True","False"],help="whether genes with 0 expression value would be included. Default 'False' is to include them")
	group3.add_argument("-cor","--correlation",default="pearson",choices=[False,"pearson","spearman"],help="select the type of correlation in the table, default is pearson")
	group2 = parser.add_argument_group('Chart visulaization arguments')
	group2.add_argument("-mean","--showmeans",default="True",choices=["True","False"],help="whether to show the position of mean in boxplot or violin plot, default 'True' is to show")
	group2.add_argument("-sf","--showfliers",default="True",choices=["True","False"],help="whether to show outliers in boxplots, default 'True' is to show")
	group2.add_argument("-ylim","--ylimit",default="False",choices=["True","False"], help="whether to show the y-axis to 100, default False is automatically adjusted")
	group4 = parser.add_argument_group('Graphing arguments')
	group4.add_argument("--dotsize",default=20,type=int,help="dotsize, default is 20")
	group4.add_argument("--textsize",default=20,type=int,help="textsize, default is 20")
	group4.add_argument("--ticksize",default=15,type=int,help="ticksize, default is 15")
	group4.add_argument("--labelsize",default=20,type=int,help="labelsize, default is 20")
	group4.add_argument("--titlesize",default=20,type=int,help="titlesize, default is 20")
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


def groupingplot(All_value,name,plot="boxplot",context="CG",target="Gene_Body",threshold=1000,skip0=False,n=3,showmeans=True,showfliers=True,ylim=False,table=True,corr="pearson",figuresize=(8,6),dotsize=30,textsize=20,ticksize=15,labelsize=20,titlesize=20):
	con=[]
	con.append(context)
	tar=[]
	tar.append(target)
	order=["1st","2nd","3rd"]
	for i in range(4,n+1):
		order.append("%dth"%(i))
	for i in con:	##["CG","CHG","CHH"]
		for j in tar:  ##["Gene_Body","Promoter","Exon","Intron"]
			groupname=order[:n]	
			sepvalue=[0]  
			groups={}  
			no_exp=All_value["%s"%(i)].loc[All_value["%s"%(i)]["RPKM"]==0]  
			other_exp=All_value["%s"%(i)].loc[All_value["%s"%(i)]["RPKM"]!=0]  
			sort_other_exp=other_exp.sort_values(by='RPKM')
			dataforplot=[]
			maxall=[]
			numerical_equation=""
			plt.figure(figsize=(8,6))
			plt.subplot(111)
			if skip0==False:
				groupname.insert(0,"no") 
				groups["no"]=no_exp
				no_exp_meth=groups["no"]["%s"%(j)].dropna() 
				dataforplot.append(no_exp_meth) 
				numerical_equation+="no (RPKM=0.00)"  
			for k in range(n):
				sepvalue.append(sort_other_exp.iloc[(len(sort_other_exp)*(k+1)//n)-1,5]) 
				groups[order[k]]=sort_other_exp.loc[(sort_other_exp["RPKM"]>sepvalue[k]) & (sort_other_exp["RPKM"]<=sepvalue[k+1])]
				yes_exp_meth=groups[order[k]]["%s"%(j)].copy().dropna()
				dataforplot.append(yes_exp_meth)
			if ylim == "True":  
				plt.ylim((0, 100))
			if plot == "boxplot":
				bp=plt.boxplot(dataforplot,patch_artist=True,showmeans=showmeans,showfliers=showfliers) 
				plt.setp(bp["means"],marker="o",markeredgecolor="red",markerfacecolor="red",markersize=3)
				plt.setp(bp["medians"],color="Navy")
				plt.setp(bp["whiskers"],color="Navy")
				plt.setp(bp["caps"],color="Navy")
				plt.setp(bp["fliers"],markeredgecolor="Navy")
				for box in bp["boxes"]:
					box.set(color="Navy")
					box.set(facecolor="SkyBlue")
			if plot == "violinplot":
				plt.violinplot(dataforplot,showmeans=showmeans)
			if skip0==False: 
				plt.xticks(list(np.arange(start=1,stop=n+2)),groupname,fontsize=ticksize) 
			else:
				plt.xticks(list(np.arange(start=1,stop=n+1)),groupname,fontsize=ticksize)
			plt.title("%s"%(i),fontsize=titlesize)
			plt.ylabel("%s Methylation "%(j)+"(%)",fontsize=labelsize)
			plt.yticks(fontsize=ticksize)
			plt.savefig('%s_%s_%s_%s_%dexpress.png'%(name,i,j,plot,n), dpi=300)
			#plt.show()	
			##########   table for grouping   #############
			corr="pearson"
			table=pd.DataFrame()
			TableList=[]
			if skip0==False:
				sepvalue.insert(0,0)  
			for k in range(len(groupname)):
				info=dataforplot[k].describe()
				RPKM_range="%-3.2f<RPKM<=%-3.2f"%(sepvalue[k],sepvalue[k+1])
				corrdata=groups[groupname[k]][["%s"%(j),"RPKM"]]
				corrdata=corrdata.dropna(axis=0,how="any")
				ar=corrdata.values
				if corr=="pearson" or "spearman" or "leastSquares":
					if corr == "pearson":
						rho, pval = scipy.stats.pearsonr(ar[:,0], ar[:,1])
					if corr == "spearman":
						rho, pval = scipy.stats.spearmanr(ar,axis=0) 
					if corr == "leastSquares":
						slope, intercept, rho, pval, std_err = scipy.stats.C(ar[:,0], ar[:,1])
				lst=[groupname[k],len(groups[groupname[k]]),int(info[0]),RPKM_range,rho,pval,info[1],info[2],info[3],info[4],info[5],info[6],info[7]]
				TableList.append(lst)
			table=table.append(TableList,ignore_index=True)
			table.columns=[["groupname","n_genes","n_dpna_genes","RPKM_range","correlation","p_value","mean","std","min","Q1","median","Q3","max"]]
			table.to_csv("%s_%s_%s_expression_%dgroup_summary.txt"%(name,i,j,n),header=True,index=False,sep="\t")
			print("complete %s on %s"%(context,target))


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
            groupingplot(All_value,name, plot=args.plot, context=con, target=tar, skip0=eval(args.skip0), n=args.numberofgroup,showmeans=eval(args.showmeans),showfliers=eval(args.showfliers),ylim=args.ylimit,corr=args.correlation,figuresize=(8,6),dotsize=args.dotsize,textsize=args.textsize,ticksize=args.ticksize,labelsize=args.labelsize,titlesize=args.titlesize)



if __name__ == '__main__':
    main()

