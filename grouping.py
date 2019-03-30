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
	group1.add_argument("-p","--plot",type=str,default="boxplot",choices=["boxplot","violinplot"],help="Create boxplot or violinplot, default is boxplot")
	group1.add_argument("-c","--context",type=str,default="CG",choices=["CG","CHG","CHH"],help="Choose the context of methylation")
	group1.add_argument("-t","--target",type=str,default="Promoter",choices=["Gene_Body","Promoter","Exon","Intron"],help="Choose the target region of methylation")
	group1.add_argument("-nb","--numberofgroup",default=5,type=int,help="Define how many group to seperate gene expression, default is 5")
	group3 = parser.add_argument_group('Important general arguments')
	#group3.add_argument("-thrs","--threshold",default="1000",help="Whether to skip genes with expression value that is too high, default is to skip genes higher than 1000 expression")
	group3.add_argument("-re0","--skip0",default="False",choices=["True","False"],help="Whether genes with 0 expression value would be included. Default is to include them")
	group3.add_argument("-cor","--correlation",default="pearson",choices=[False,"pearson","spearman"],help="select the type of correlation in the table")
	group2 = parser.add_argument_group('Chart visulaization arguments')
	group2.add_argument("-mean","--showmeans",default="True",choices=["True","False"],help="whether to show the position of mean in boxplot or violin plot, default is to show")
	group2.add_argument("-sf","--showfliers",default="True",choices=["True","False"],help="remove the outliers in the boxplots")
	group2.add_argument("-ylim","--ylimit",default="False",help="Nemeric value. Optional. zoom in the DNA methylation level to clearly understand the methylation distribution") ##可以改
	#parser.add_argument("-table","--createtable",default="True",choices=["True","False"],help="create table to show the statistics in each group")##False是否可以去掉
	group4 = parser.add_argument_group('Graphing arguments')
	#parser.add_argument("-fig","--figuresize",default=(8,6),help="the figure size")
	group4.add_argument("--dotsize",default=20,type=int,help="dotsize")
	group4.add_argument("--textsize",default=20,type=int,help="textsize")
	group4.add_argument("--ticksize",default=15,type=int,help="ticksize")
	group4.add_argument("--labelsize",default=20,type=int,help="labelsize")
	group4.add_argument("--titlesize",default=20,type=int,help="titlesize")
	group4.add_argument("--legendsize",default=20,type=int,help="legendsize")
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
			groupname=order[:n]	## for xticks usage, should put in the loop
			sepvalue=[0]  ##區分exp value的那個值，第一個值是0
			groups={}  ##用來放不同群gene的dataframe
			no_exp=All_value["%s"%(i)].loc[All_value["%s"%(i)]["RPKM"]==0]  ##value=0
			other_exp=All_value["%s"%(i)].loc[All_value["%s"%(i)]["RPKM"]!=0]  ##vlaue!=0
			sort_other_exp=other_exp.sort_values(by='RPKM')
			##################   threshold should add in here，另一種模式的
			dataforplot=[]
			maxall=[]
			numerical_equation=""
			plt.figure(figsize=(8,6))
			plt.subplot(111)
			if skip0==False:
				groupname.insert(0,"no")  ## 分群的名稱，要加no expression
				groups["no"]=no_exp
				no_exp_meth=groups["no"]["%s"%(j)].dropna() 
				dataforplot.append(no_exp_meth) ## the data for plotting is a list but 
				#maxall.append(max(no_exp_meth))  ## put the each group value max value in the list
				numerical_equation+="no (RPKM=0.00)"  ##should become a dataframe
			for k in range(n):
				sepvalue.append(sort_other_exp.iloc[(len(sort_other_exp)*(k+1)//n)-1,5]) ###除法是怎麼算的，>>好像是無條件捨去？？？？　決定求出位置後減一了,5是gene expression value
				groups[order[k]]=sort_other_exp.loc[(sort_other_exp["RPKM"]>sepvalue[k]) & (sort_other_exp["RPKM"]<=sepvalue[k+1])]
				yes_exp_meth=groups[order[k]]["%s"%(j)].copy().dropna()  ##
				#maxall.append(max(yes_exp_meth))
				dataforplot.append(yes_exp_meth)
				#numerical_equation+=", %s (%-3.2f<RPKM<=%-3.2f)"%(order[:n][k],sepvalue[k],sepvalue[k+1])
			#highest_value=max(maxall)
			#ylimit=highest_value*3/2
			if ylim == "True":  ##另一種是此>>>
				plt.ylim((0, 100))
			if plot == "boxplot":
				bp=plt.boxplot(dataforplot,patch_artist=True,showmeans=showmeans,showfliers=showfliers) ##直接用removout ##showmeans=True>>this parameter shows the qualtiles and the means #patch_artist=True: A patch is a 2D artist with a face color and an edge color.
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
				#ax = sns.violinplot(data=dataforplot)
			if skip0==False: ##可是showno用了兩次if 
				plt.xticks(list(np.arange(start=1,stop=n+2)),groupname,fontsize=ticksize) ##list(np.arange(start=1,stop=n+2)) is like [1,2,3,4,5],,
			else:
				plt.xticks(list(np.arange(start=1,stop=n+1)),groupname,fontsize=ticksize)
			plt.title("%s"%(i),fontsize=titlesize)
			plt.ylabel("%s Methylation "%(j)+"(%)",fontsize=labelsize)
			plt.yticks(fontsize=ticksize)
			#plt.xlabel("%s"%(numerical_equation),fontsize=9)
			plt.savefig('%s_%s_%s_%s_%dexpress.png'%(name,i,j,plot,n))
			plt.show()	
			##########   table for grouping   #############
			corr="pearson"
			table=pd.DataFrame()##,columns=["groupname","n_genes","n_nan_genes","RPKM_range","correlatioin","p_value","max","min","mean","std"],dtype={"groupname":str,"n_genes":int,"RPKM_range":str,"correlatioin":float,"max":float,"min":float,"mean":float,"sd":float}
			#########   table for plotting
			TableList=[]
			if skip0==False:
				sepvalue.insert(0,0)  ##################################################  這放在迴圈裡很危險,但若放在上面的迴圈就還好
			for k in range(len(groupname)):
				info=dataforplot[k].describe()
				RPKM_range="%-3.2f<RPKM<=%-3.2f"%(sepvalue[k],sepvalue[k+1])
				corrdata=groups[groupname[k]][["%s"%(j),"RPKM"]]
				corrdata=corrdata.dropna(axis=0,how="any")
				ar=corrdata.values
				if corr=="pearson" or "spearman" or "leastSquares":
					if corr == "pearson":
						rho, pval = scipy.stats.pearsonr(ar[:,0], ar[:,1]) ##pearson
						#plt.text(1,1,"R = %2.2f\nP = %2.2f"%(rho,pval),fontsize=10,horizontalalignment='right',verticalalignment='top')
					if corr == "spearman":
						rho, pval = scipy.stats.spearmanr(ar,axis=0) ##spearman , axis=0
						#plt.text(max(ar[:,0]),max(ar[:,1]),"R = %2.2f\nP = %2.2f"%(rho,pval),fontsize=10,horizontalalignment='right',verticalalignment='top')
					if corr == "leastSquares":
						slope, intercept, rho, pval, std_err = scipy.stats.C(ar[:,0], ar[:,1])
				lst=[groupname[k],len(groups[groupname[k]]),int(info[0]),RPKM_range,rho,pval,info[1],info[2],info[3],info[4],info[5],info[6],info[7]]
				TableList.append(lst)
			table=table.append(TableList,ignore_index=True)
			table.columns=[["groupname","n_genes","n_dpna_genes","RPKM_range","correlation","p_value","mean","std","min","Q1","median","Q3","max"]]
			print table
			table.to_csv("%s_%s_%s_expression_%dgroup_summary.txt"%(name,i,j,n),header=True,index=False,sep="\t")


def main():
    parser = get_parser()
    args = parser.parse_args()
    All_value,name=read_meth_exp_file(args.samplename)
    groupingplot(All_value,name,plot=args.plot,context=args.context,target=args.target,skip0=eval(args.skip0),n=args.numberofgroup,showmeans=eval(args.showmeans),showfliers=eval(args.showfliers),ylim=eval(args.ylimit),corr=args.correlation,figuresize=(8,6),dotsize=args.dotsize,textsize=args.textsize,ticksize=args.ticksize,labelsize=args.labelsize,titlesize=args.titlesize)#table=eval(args.createtable),,threshold=eval(args.threshold)



if __name__ == '__main__':
    main()

