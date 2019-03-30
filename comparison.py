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
    parser.add_argument("-s","--samplelist",default = "samplelist.txt",type=str,help="put in the sample description file")
    parser.add_argument("-p","--plot",type=str,default="scatter",choices=["scatter","heatmap"],help="create scatterplot or heatmap")
    parser.add_argument("-c","--context",type=str,default="CG",choices=["CG","CHG","CHH"],help="choose the context of methylation")
    parser.add_argument("-t","--target",type=str,default="Gene_Body",choices=["Gene_Body","Promoter","Exon","Intron"],help="choose the target region of methylation")
    parser.add_argument("-mthr","--meththreshold",default=10,type=float,help="set threshold to identify the differential methylated genes")
    parser.add_argument("-ethr","--expthreshold",default=1,type=float,help="set threshold to identify genes that have expression change")
    group3 = parser.add_argument_group('Important general arguments')
    ##parser.add_argument("-thrs","--threshold",default="1000",help="set the reasonable value of gene expression value to avoid outliers")####>>>>>>>>>>>>???
    #parser.add_argument("-re0","--skip0",default="False",choices=["True","False"],help="remove genes that their expression value is equal to 0") ####>>>>>>???
    parser.add_argument("-ad","--addGEvalue",default=1.0,type=float,help="add a small value on gene expression value to calculate the log(fold change)")
    ##parser.add_argument("-ylim","--ylimit",default="False",help="zoom in the DNA methylation value to clearly understand the distribution")
    parser.add_argument("-list","--genelist",default="True",choices=["True","False"],help="create table to show information of the genes selected")
    parser.add_argument("-cor","--correlation",default="pearson",choices=["False","pearson","spearman"],help="select the type of correlation, default is 'pearson'")
    group4 = parser.add_argument_group('Graphing arguments')
    #parser.add_argument("-fig","--figuresize",default=(8,6),help="the figure size")
    group4.add_argument("--dotsize",default=20,type=int,help="dotsize")
    group4.add_argument("--textsize",default=25,type=int,help="textsize")
    group4.add_argument("--ticksize",default=15,type=int,help="ticksize")
    group4.add_argument("--labelsize",default=25,type=int,help="labelsize")
    group4.add_argument("--titlesize",default=25,type=int,help="titlesize")
    group4.add_argument("--legendsize",default=20,type=int,help="legendsize")
    group4.add_argument("--fontsize",default=1.2,type=int,help="fontsize")
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
    return All_value




def compare_value(Samplelist="samplelist.txt",context="CG",target="Gene_Body",add=1,Methcutoff=10,Expcutoff= 1, savelist=True):
	Samplelist=pd.read_csv("%s"%(Samplelist),header=None,sep="\t") ##names,cgmap,exp,trt
	##這三行就把所有直放進去了
	con=[]
	con.append(context)
	tar=[]
	tar.append(target)
	value={}
	for i in con: ##["CG","CHG","CHH"]
		for j in tar: ##["Gene_Body","Promoter","Exon","Intron"]
			for k in Samplelist[0]:  ## 用k 就好了吧
				value[k]=read_meth_exp_file(k)  ## 裡面是sample name
			meanvalue=value[k][i][["gene_ID"]]  ##先放入ｇｅｎｅ＿ＩＤ  ##because of k, so should contain the k in value
			everyvalue=value[k][i][["gene_ID"]]
			##methexpall>>>add 4meth, 4exp  , just add formethaverage, forexpaverage in the dataframe
			treatment=Samplelist[3].unique()
			##create meanvalue for comparison
			for trt in treatment:  ###list different treatment name 
				formethaverage=pd.DataFrame() ## empty
				forexpaverage=pd.DataFrame() ## empty
				Samples=Samplelist.loc[Samplelist[3]==trt] ## the data from only that treatment
				for k in range(len(Samples)):
					formethaverage[Samples.iloc[k,0]]=value[Samples.iloc[k,0]][i][j] 
					everyvalue["%smeth"%(Samples.iloc[k,0])]=value[Samples.iloc[k,0]][i][j]
					forexpaverage[Samples.iloc[k,0]]=value[Samples.iloc[k,0]][i]["RPKM"]
					everyvalue["%sexp"%(Samples.iloc[k,0])]=value[Samples.iloc[k,0]][i]["RPKM"]
				meanvalue["%smethmean"%(trt)]=formethaverage.mean(1).copy()  ## average between samples  
				meanvalue["%sexpmean"%(trt)]=forexpaverage.mean(1).copy()+add ## modify RPKM value to avoid devide 0 error
			dropres=meanvalue.dropna(axis=0,how="any")  ###drop the value that is zero ##33357-33020
			print "genes that have values in both treatment", len(dropres)
			####  minus the average value (後減前MT-WT or cancer-normal)
			dropres["Δmethylation"]=dropres["%smethmean"%(treatment[1])].copy()-dropres["%smethmean"%(treatment[0])].copy()
			dropres["Δgene expression"]=dropres["%sexpmean"%(treatment[1])].copy()/dropres["%sexpmean"%(treatment[0])].copy()  ##相除
			dropres["Δgene expression"]=dropres["Δgene expression"].apply(np.log2)
			## start to cut off 
			selectMeth=dropres.loc[(abs(dropres["Δmethylation"])>Methcutoff)]  ##31>>>can use to draw heatmap
			print "genes that are differential methylated between treatment",len(selectMeth)
			selectMeth.reset_index(inplace=True) ## to let the index in order
			##select meth&exp
			selectMethExp=selectMeth.loc[(abs(selectMeth["Δgene expression"])>Expcutoff)]  ##31
			print "genes that are differential expressed",len(selectMethExp)
			if savelist:
				### save genes
				selectMethExp.iloc[:,1:].to_csv("chosen_genes_meth%.2f_exp%.2f.txt"%(Methcutoff,Expcutoff),sep="\t",index=False,header=True)
	return Samplelist,dropres,selectMeth,selectMethExp



def CompareScatter(dropres,selectMethExp,context="CG",target="Gene_Body",corr="pearson",figuresize=(8,6.5),dotsize=20,textsize=20,ticksize=15,labelsize=20,titlesize=20):
	con=[]
	con.append(context)
	tar=[]
	tar.append(target)
	value={}
	for i in con: ##["CG","CHG","CHH"]
		for j in tar: ##["Gene_Body","Promoter","Exon","Intron"]	
			fig=plt.figure(figsize=figuresize)
			ax1=fig.add_subplot(111)
			ax1.scatter(list(dropres["Δmethylation"]),list(dropres["Δgene expression"]),s=dotsize) ##size is s
			##add select plots
			ax1.scatter(list(selectMethExp["Δmethylation"]),list(selectMethExp["Δgene expression"]),s=dotsize,c="red")  ##red point cover the old one
			ax1.set_title("%s %s"%(i,j),fontsize=titlesize)
			while corr=="pearson" or "spearman" or "leastSquares":
				if corr == "pearson":
					rho, pval = scipy.stats.pearsonr(dropres["Δmethylation"], dropres["Δgene expression"]) ##pearson
					#plt.text(1,1,"R = %2.2f\nP = %2.2f"%(rho,pval),fontsize=10,horizontalalignment='right',verticalalignment='top')
				if corr == "spearman":
					rho, pval = scipy.stats.spearmanr(ar,axis=0) ##spearman , axis=0
					#plt.text(max(ar[:,0]),max(ar[:,1]),"R = %2.2f\nP = %2.2f"%(rho,pval),fontsize=10,horizontalalignment='right',verticalalignment='top')
				if corr == "leastSquares":
					slope, intercept, rho, pval, std_err = scipy.stats.C(ar[:,0], ar[:,1])
				a=max(list(dropres["Δmethylation"]))
				b=max(list(dropres["Δgene expression"]))
				plt.text(a,b,"R = %-2.3f\nP = %-2.2e"%(rho,pval),fontsize=textsize,horizontalalignment='right',verticalalignment='top') ##max(ar[:,0]),ylimit*0.9,
			plt.xlabel(r'$\Delta\ methylation$'+" (%)",fontsize=labelsize)
			plt.ylabel(r'$\Delta\ gene\ expression\ (log2\ fold\ change)$',fontsize=labelsize) ##need blank"\ "
			plt.xticks(fontsize=ticksize)
			plt.yticks(fontsize=ticksize)
			plt.axhline(0,color="gray",linestyle="--") ##水平
			#plt.axhline(1,color="peru",linestyle="--",linewidth=0.5)
			#plt.axhline(-1,color="peru",linestyle="--",linewidth=0.5)
			plt.axvline(0,color="gray",linestyle="--") ##垂直
			#plt.axvline(20,color="peru",linestyle="--",linewidth=0.5) 
			#plt.axvline(-20,color="peru",linestyle="--",linewidth=0.5) 
			plt.savefig('%s_%s_comparisonplot.png'%(i,j))
			#plt.plot([-20,],[])
			plt.show()



def heatmap(Samplelist,selectMeth,context="CG",target="Gene_Body",maxexp=20,Expcutoff=2,fontsize=1.2):
	con=[]
	con.append(context)
	tar=[]
	tar.append(target)
	value={}
	for i in con: ##["CG","CHG","CHH"]
		for j in tar: ##["Gene_Body","Promoter","Exon","Intron"]
			treatment=Samplelist[3].unique()
			##從selectMeth開始>>>將 exp大於多少的都算是那個值
			selectMeth["WTexpmean"].loc[selectMeth["WTexpmean"]>maxexp]=maxexp
			selectMeth["OTU5expmean"].loc[selectMeth["OTU5expmean"]>maxexp]=maxexp
			###  to let the exp value from (0,20) to (0,100)
			selectMeth["WTexpmean"]=selectMeth["WTexpmean"]*5
			selectMeth["OTU5expmean"]=selectMeth["OTU5expmean"]*5
			##這有warning喔
			### trancate first than decide the cutoff of differential expression
			selectMeth["Δgene expression"]=selectMeth["%sexpmean"%(treatment[1])].copy()/selectMeth["%sexpmean"%(treatment[0])].copy()  ##相除
			selectMethExp=selectMeth.loc[(abs(selectMeth["Δgene expression"])>Expcutoff)]
			print "selectmeth",len(selectMeth)
			print "genes that are differential expressed",len(selectMethExp)
			#p=sns.clustermap(selectMeth[["WTmethmean","OTU5methmean"]])##,"WTexp","OTU5exp">>>>改成只用cluster就好
			sns.set(font_scale=fontsize)## scale all fonts in your legend and on the axes.
			##改成meth,exp 有差異都畫
			p=sns.clustermap(selectMethExp[["WTmethmean","OTU5methmean","WTexpmean","OTU5expmean"]],vmin=0,vmax=100,robust=True,col_cluster=False,yticklabels=False)  ##這樣clustering 比較漂亮  ##,center=10
			##p.ax_heatmap.set_title("%s %s"%(i,j),fontsize=titlesize)##plt.title("%s %s"%(i,j),fontsize=titlesize)>>add to legend
			plt.savefig('%s_%s_heatmap.png'%(i,j))
			plt.show()
	return selectMethExp



def main():
    parser = get_parser()
    args = parser.parse_args()
    Samplelist,dropres,selectMeth,selectMethExp=compare_value(Samplelist=args.samplelist,context=args.context,target=args.target,add=args.addGEvalue,Methcutoff=args.meththreshold,Expcutoff= args.expthreshold)
    if args.plot=="scatter" :
    	CompareScatter(dropres,selectMethExp,corr=args.correlation,context=args.context,target=args.target,figuresize=(10,8),dotsize=args.dotsize,textsize=args.textsize,ticksize=args.ticksize,labelsize=args.labelsize,titlesize=args.titlesize)
    if args.plot=="heatmap" :
    	heatmapvalue=heatmap(Samplelist,selectMeth,context=args.context,target=args.target,maxexp=20,Expcutoff=args.expthreshold,fontsize=args.fontsize)



if __name__ == '__main__':
    main()




