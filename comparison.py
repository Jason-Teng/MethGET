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
import time
import os, sys
import scipy
import seaborn as sns
import warnings
from sklearn.mixture import GaussianMixture
from scipy.stats import norm, multivariate_normal
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
    group1.add_argument("-t","--target",type=str,default="Promoter",choices=["Promoter","Gene_Body","Exon","Intron","all"],help="choose the genomic location of methylation, default is Promoter")
    group1.add_argument("-pro","--prob_cutoff", default="0.000001", type=float, help="define the differential genes by their probrability in Gaussian mixture model, default is '0.000001'")
    group2 = parser.add_argument_group('Important general arguments')
    group2.add_argument("-p","--plot",type=str,default="scatter",choices=["scatter","kernel"],help="create scatterplot or kernel density plot, default is scatter")
    group2.add_argument("-cor","--correlation",default="pearson",choices=["False","pearson","spearman"],help="select the type of correlation, default is pearson")
    group2.add_argument("--shownumber",default="False",choices=["False","True"],help="whether to show the number of significant genes, default False is not show")
    group3 = parser.add_argument_group("Define anomaly genes with methylation changes and gene expression changes")
    group3.add_argument("--cutoff", default="False", choices=["False","True"], help="whether to use methylation and gene expression cutoff to define outliers, default is False")
    group3.add_argument("-mthr","--meththreshold",default="auto",help="set cutoff of differential methylated genes. default 'auto' uses methylation changes, CG:10, CHG:1, CHH:1")
    group3.add_argument("-ethr","--expthreshold",default=1.0,type=float,help="set cutoff to identify genes that have expression changes, default uses expression changes: log2FC is 1")
    group3.add_argument("--methmin",default=None,help="minimum methylation changes for x-axis")
    group3.add_argument("--methmax",default=None,help="maximum methylation changes for x-axis")
    group3.add_argument("--expmin",default=None,help="minimum gene expression changes for y-axis")
    group3.add_argument("--expmax",default=None,help="maximum gene expression changes for y-axis")
    group4 = parser.add_argument_group('Graphing arguments')
    group4.add_argument("--dotsize",default=20,type=int,help="dotsize, default is 20")
    group4.add_argument("--textsize",default=25,type=int,help="textsize, default is 25")
    group4.add_argument("--ticksize",default=15,type=int,help="ticksize, default is 15")
    group4.add_argument("--labelsize",default=25,type=int,help="labelsize, default is 25")
    group4.add_argument("--titlesize",default=25,type=int,help="titlesize, default is 25")
    group4.add_argument("--fontsize",default=1.2,type=float,help="fontsize, default is 1.2")
    return parser



def read_meth_exp_file(name):
    """
    read all the preprocess files by the sample name
    """
    exp=pd.read_csv("%s_exp.txt"%(name),sep="\t",names=["gene_ID","RPKM"]) # read expression file
    GeneMeth=pd.read_csv("%s_GeneMean.txt"%(name),sep="\t",names=["gene_ID","Meth_value_mean_CG","Meth_value_mean_CHG","Meth_value_mean_CHH"]) # read gene body methylation file
    PromoterMeth=pd.read_csv("%s_PromoterMean.txt"%(name),sep="\t",names=["gene_ID","Meth_value_mean_CG","Meth_value_mean_CHG","Meth_value_mean_CHH"]) # # read promoter methylation file
    ExonMeth=pd.read_csv("%s_ExonMean.txt"%(name),sep="\t",names=["gene_ID","Meth_value_mean_CG","Meth_value_mean_CHG","Meth_value_mean_CHH"]) # read exon methylation file
    IntronMeth=pd.read_csv("%s_IntronMean.txt"%(name),sep="\t",names=["gene_ID","Meth_value_mean_CG","Meth_value_mean_CHG","Meth_value_mean_CHH"]) # read intron methylation file
    # to combine all files into a dict (contexts and genomic locations)
    All_value={}
    for i in ["CG","CHG","CHH"]:
        meth=pd.DataFrame()
        meth["gene_ID"]=GeneMeth["gene_ID"]
        meth["Gene_Body"]=GeneMeth["Meth_value_mean_%s"%(i)]*100
        meth["Promoter"]=PromoterMeth["Meth_value_mean_%s"%(i)]*100
        meth["Exon"]=ExonMeth["Meth_value_mean_%s"%(i)]*100
        # merge the intron files because maybe there are genes without intron
        methin=pd.merge(meth,IntronMeth[["gene_ID","Meth_value_mean_%s"%(i)]],on=["gene_ID"],how="outer")
        intronname="Meth_value_mean_%s"%(i)
        methin=methin.rename(columns={intronname: 'Intron'})
        methin["Intron"]=methin["Intron"]*100
        # merge expression files (only genes with expression values)
        res=pd.merge(methin,exp,on=["gene_ID"],how="inner") 
        All_value["%s"%(i)]=res
    return All_value,name




def compare_value(Samplelist="samplelist.txt", context="CG", target="Gene_Body", add=1, prob_cutoff = 0.000001, usecutoff="False", Methcutoff="auto", Expcutoff=2):
    """
    Get the differential genes
    """
    Samplelist=pd.read_csv(Samplelist,header=None,sep="\t") # read samplelist file
    con=[]
    con.append(context) # context
    tar=[]
    tar.append(target) # genomic location
    value={}
    for i in con: 
        for j in tar:
            # retrieve data for certain context and genomic location
            for k in Samplelist[0]: 
                value[k],name=read_meth_exp_file(k)  
            meanvalue=value[k][i][["gene_ID"]]  
            everyvalue=value[k][i][["gene_ID"]]
            treatment=Samplelist[3].unique()
            ## create meanvalue for comparison
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
            dropres=meanvalue.dropna(axis=0,how="any") # drop genes with nan value
            dropres["Δmethylation"]= dropres["%smethmean"%(treatment[0])].copy()-dropres["%smethmean"%(treatment[1])].copy() #methylation changes
            dropres["Δgene expression"]=dropres["%sexpmean"%(treatment[0])].copy()/dropres["%sexpmean"%(treatment[1])].copy() #expression changes
            dropres["Δgene expression"]=dropres["Δgene expression"].apply(np.log2) # log2 fold change
            ## calulate the variation
            deltaMethmean = dropres["Δmethylation"].mean()
            deltaExpmean = dropres["Δgene expression"].mean()
            dropres["var"] = (dropres["Δmethylation"]-deltaMethmean)**2 + (dropres["Δgene expression"]-deltaExpmean)**2
            variation = dropres["var"].sum()/(len(dropres)-2)
            sqrtvar = np.sqrt(variation)
            ## using cutoff to define differential genes
            if usecutoff =="True":
                # the default cutoff of contexts
                if Methcutoff == "auto":
                    if i == "CG":
                        Methcutoff = float(10)
                    if i =="CHG":
                        Methcutoff = float(1)
                    if i == "CHH":
                        Methcutoff = float(1) 
                else:
                    Methcutoff = float(Methcutoff)
                log2Expcutoff = Expcutoff # the user need to provide lo2 FC for the cutoff of differential genes
                selectMethExp_1 = dropres.loc[(dropres["Δmethylation"]>Methcutoff) & (dropres["Δgene expression"]>log2Expcutoff)] # genes in 1st quardrant
                selectMethExp_2 = dropres.loc[(dropres["Δmethylation"]<Methcutoff*(-1)) & (dropres["Δgene expression"]>log2Expcutoff)] # genes in 2nd quardrant
                selectMethExp_3 = dropres.loc[(dropres["Δmethylation"]<Methcutoff*(-1)) & (dropres["Δgene expression"]<log2Expcutoff*(-1))] # genes in 3rd quardrant
                selectMethExp_4 = dropres.loc[(dropres["Δmethylation"]>Methcutoff) & (dropres["Δgene expression"]<log2Expcutoff*(-1))] # genes in 4th quardrant
                selectMeth=dropres.loc[(abs(dropres["Δmethylation"])>Methcutoff)] # genes with methylation changes
                selectMeth.reset_index(inplace=True)  # reindex the row
                selectMethExp=selectMeth.loc[(abs(selectMeth["Δgene expression"])>log2Expcutoff)]  # genes with methylation and expression changes
            ## using GMM to define differential
            if usecutoff=="False":                
                selectMeth = []
                X = np.array(dropres[["Δmethylation", "Δgene expression"]]) # data for plotting 
                X_size = X.shape[0]
                groups = 1
                gmm = GaussianMixture(n_components=groups,covariance_type='full').fit(X) # GMM model
                b=gmm.score_samples(X) # average log-likelihood
                dropres["GMM_score"] =b 
                probrability_cutoff = prob_cutoff # the threshold of likelihood
                log_prob_cutoff = np.log(prob_cutoff) # log likelihood
                selectMethExp = dropres.loc[(dropres["GMM_score"]<log_prob_cutoff)] # to select differential genes with low likelihood
                selectMethExp_1 = selectMethExp[(dropres["Δmethylation"]>0) & (dropres["Δgene expression"]>0)] # genes in 1st quardrant
                selectMethExp_2 = selectMethExp[(dropres["Δmethylation"]<0) & (dropres["Δgene expression"]>0)] # genes in 2nd quardrant
                selectMethExp_3 = selectMethExp[(dropres["Δmethylation"]<0) & (dropres["Δgene expression"]<0)] # genes in 3rd quardrant
                selectMethExp_4 = selectMethExp[(dropres["Δmethylation"]>0) & (dropres["Δgene expression"]<0)] # genes in 4th quardrant
    return Samplelist,dropres,selectMeth,selectMethExp, selectMethExp_1, selectMethExp_2, selectMethExp_3, selectMethExp_4, deltaMethmean, deltaExpmean, variation


def CompareScatterplot(dropres,selectMethExp,selectMethExp_1,selectMethExp_2,selectMethExp_3,selectMethExp_4, variation, plot="scatter",context="CG",target="Gene_Body",corr="pearson",shownumber=True,methmin=None,methmax=None,expmin=None,expmax=None,figuresize=(10,8),dotsize=15,textsize=20,ticksize=15,labelsize=20,titlesize=20):
    """
    generate the scatterplot for comparison analyses
    """
    con=[]
    con.append(context) # contexts
    tar=[]
    tar.append(target) # genomic location
    value={}
    for i in con: ##["CG","CHG","CHH"]
        for j in tar: ##["Gene_Body","Promoter","Exon","Intron"]    
            ### scatter plot
            if plot == "scatter":
                fig=plt.figure(figsize=figuresize)
                ax1=fig.add_subplot(111)
                ax1.scatter(list(dropres["Δmethylation"]),list(dropres["Δgene expression"]),s=dotsize) # plot all genes
                ax1.scatter(list(selectMethExp["Δmethylation"]),list(selectMethExp["Δgene expression"]),s=dotsize,c="red") # plot differential genes           
            # the position to set the text of the figure
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
            ### kernel density plot
            if plot == "kernel":
                plt.figure(figsize=figuresize)
                ax1=sns.kdeplot(list(dropres["Δmethylation"]),list(dropres["Δgene expression"]),cmap="Blues", shade=True)
            ## the correlation coefficient
            while corr=="pearson" or corr=="spearman" or corr=="leastSquares":
                if corr == "pearson":
                    rho, pval = scipy.stats.pearsonr(dropres["Δmethylation"], dropres["Δgene expression"])
                if corr == "spearman":
                    rho, pval = scipy.stats.spearmanr(dropres["Δmethylation"], dropres["Δgene expression"])
                if corr == "leastSquares":
                    slope, intercept, rho, pval, std_err = scipy.stats.C(ar[:,0], ar[:,1])
                plt.text(a,b,"R = %-2.3f\nP = %-2.2e"%(rho,pval),fontsize=textsize,horizontalalignment='right',verticalalignment='top') # plot the correlation coefficient
                break
            ax1.set_title("%s %s"%(i,j),fontsize=titlesize)
            if shownumber == "True":
                ## show the differential genes in each quardrants
                plt.text(a,b,"%d"%(len(selectMethExp_1)),fontsize=textsize,horizontalalignment='right',verticalalignment='top')
                plt.text(c,b,"%d"%(len(selectMethExp_2)),fontsize=textsize,horizontalalignment='left',verticalalignment='top')
                plt.text(c,d,"%d"%(len(selectMethExp_3)),fontsize=textsize,horizontalalignment='left',verticalalignment='bottom')
                plt.text(a,d,"%d"%(len(selectMethExp_4)),fontsize=textsize,horizontalalignment='right',verticalalignment='bottom')
            ## graphing parameters 
            plt.xlabel(r'$Methylation\ level\ changes$' + " (%)",fontsize=labelsize,)
            plt.ylabel(r'$Gene\ expression\ changes\ (log2\ FC)$',fontsize=labelsize)
            plt.xticks(fontsize=ticksize)
            plt.yticks(fontsize=ticksize)
            plt.xlim(methmin,methmax)
            plt.ylim(expmin,expmax)
            plt.axhline(0,color="gray",linestyle="--")
            plt.axvline(0,color="gray",linestyle="--") 
            deltascatterplotname = 'comparison%s_plot.png'%(plot)
            plt.savefig(i + "_" + j +"_"+ deltascatterplotname, dpi = 300) # save the figure
            ## the information of differential genes in the table
            deltascattertablename = "chosen_genes.txt" 
            try:
                ratio = float(len(selectMethExp_1)*len(selectMethExp_3))/(len(selectMethExp_2)*len(selectMethExp_4))
                deltainfo = "#quardrant_1: %d, #quardrant_2: %d, #quardrant_3: %d, #quardrant_4: %d."%(len(selectMethExp_1), len(selectMethExp_2), len(selectMethExp_3), len(selectMethExp_4))
            except:
                ratio = "division by zero"
                deltainfo = "#quardrant_1: %d, #quardrant_2: %d, #quardrant_3: %d, #quardrant_4: %d."%(len(selectMethExp_1), len(selectMethExp_2), len(selectMethExp_3), len(selectMethExp_4))
            print("information for %s, %s"%(context, target), deltainfo)
            df = selectMethExp.iloc[:,0:-2]        
            df.to_csv(i + "_" + j + "_" + deltascattertablename, sep="\t",index=False, header=True) # save the table
            #plt.show()
            plt.close()
    return deltascatterplotname, deltascattertablename



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
    meth_con = [args.context] # context
    meth_tar = [args.target] # genomic location
    if args.context == "all":
        meth_con = ["CG", "CHG", "CHH"]
    if args.target == "all":
        meth_tar = ["Promoter","Gene_Body","Exon","Intron"]
    for con in meth_con:
        for tar in meth_tar:
            corr = args.correlation # if show correlation coefficient
            if args.shownumber==True:
                corr = False # not to show correlation on the plot
            Samplelist,dropres,selectMeth,selectMethExp, selectMethExp_1, selectMethExp_2, selectMethExp_3, selectMethExp_4, deltaMethmean, deltaExpmean, variation = compare_value(Samplelist=args.samplelist,context=con, target=tar,add=1.0, prob_cutoff = args.prob_cutoff, usecutoff=args.cutoff, Methcutoff=args.meththreshold,Expcutoff= args.expthreshold) # get the values of all genes and differential genes
            print("number of all genes: %d"%(len(dropres.dropna(axis=0,how="any")))) # print the number of all genes
            deltascatterplotname, deltascattertablename = CompareScatterplot(dropres, selectMethExp, selectMethExp_1, selectMethExp_2, selectMethExp_3, selectMethExp_4, variation, plot=args.plot, context=con,target=tar, corr=corr, shownumber=args.shownumber,methmin=args.methmin,methmax=args.methmax,expmin=args.expmin,expmax=args.expmax,figuresize=(10,8),dotsize=args.dotsize,textsize=args.textsize,ticksize=args.ticksize,labelsize=args.labelsize,titlesize=args.titlesize) # create the final figure
            print("number of selected genes: %d"%(len(selectMethExp))) # print the number of differential genes
    script_end = time.time() # end time
    print("preprocess time %s"%HowManyTime(script_start,script_end)) # total time for analyses



if __name__ == '__main__':
    main()

