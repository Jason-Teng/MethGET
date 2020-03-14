# coding=UTF-8

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
    group1.add_argument("-t","--target",type=str,default="Promoter",choices=["Promoter", "Gene_Body","Exon","Intron"],help="choose the genomic location of methylation, default is Promoter")
    group1.add_argument("-pro","--prob_cutoff", default="0.000001", type=float, help="define the differential genes by their probrability in Gaussian mixture model, default is '0.000001'")
    group2 = parser.add_argument_group('Important general arguments')
    group2.add_argument("-mmax","--mmax",default=100,type=float,help="set the max methylation value for heatmap, default is 100")
    group2.add_argument("-emax","--emax",default=20,type=float,help="set the max expression value for heatmap, default is 20")
    group3 = parser.add_argument_group('Define anomaly genes with methylation changes and gene expression changes')
    group3.add_argument("--cutoff", default="False", choices=["False","True"], help="whether to use methylation and gene expression cutoff to define outliers, default is False")
    group3.add_argument("-mthr","--meththreshold",default="auto",help="set cutoff of differential methylated genes. default 'auto' uses changes of methylation, CG:10, CHG:1, CHH:1")
    group3.add_argument("-ethr","--expthreshold",default=1.0,type=float,help="set cutoff to identify genes that have expression changes, default uses expression changes: log2FC is 1")
    group4 = parser.add_argument_group('Graphing arguments')
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
                selectMethExp = selectMethExp.drop("GMM_score", axis = 1)
                selectMethExp.reset_index(inplace=True)
    return Samplelist,dropres,selectMeth,selectMethExp, selectMethExp_1, selectMethExp_2, selectMethExp_3, selectMethExp_4, deltaMethmean, deltaExpmean, variation




def heatmap(Samplelist, dropres, selectMeth, selectMethExp, context="CG",target="Gene_Body",maxexp=20 ,maxmeth="auto", Expcutoff=2,fontsize=1.2):
    """
    generate heatmap figures
    """
    con=[]
    con.append(context) # context
    tar=[]
    tar.append(target) # genomic locations
    value={}
    for i in con:
        for j in tar:
            treatment = Samplelist[3].unique() # the name of the two group
            # the default maximum methylation level of contexts
            if maxmeth == "auto":
                if i == "CG":
                    maxmeth = float(100)
                if i =="CHG":
                    maxmeth = float(10)
                if i == "CHH":
                    maxmeth = float(10)
            if maxmeth != "auto":
                maxmeth = float(maxmeth)
            if True:
                ########## expression
                ## to truncate
                selectMethExp["%sexpmean"%(treatment[0])].loc[selectMethExp["%sexpmean"%(treatment[0])]>maxexp]=maxexp
                selectMethExp["%sexpmean"%(treatment[1])].loc[selectMethExp["%sexpmean"%(treatment[1])]>maxexp]=maxexp
                ###  to let the exp value from (0,20) to (0,100), multiply the scale
                selectMethExp["%sexpmean"%(treatment[0])]=selectMethExp["%sexpmean"%(treatment[0])]*(100.0/maxexp)
                selectMethExp["%sexpmean"%(treatment[1])]=selectMethExp["%sexpmean"%(treatment[1])]*(100.0/maxexp)
                ########## methylation
                ## to truncate
                selectMethExp["%smethmean"%(treatment[0])].loc[selectMethExp["%smethmean"%(treatment[0])]>maxmeth]=maxmeth
                selectMethExp["%smethmean"%(treatment[1])].loc[selectMethExp["%smethmean"%(treatment[1])]>maxmeth]=maxmeth
                ###  to let the meth value from (0,20) to (0,100), multiply the scale
                selectMethExp["%smethmean"%(treatment[0])]=selectMethExp["%smethmean"%(treatment[0])]*(100.0/maxmeth)
                selectMethExp["%smethmean"%(treatment[1])]=selectMethExp["%smethmean"%(treatment[1])]*(100.0/maxmeth)
            ### plot heatmap
            sns.set(font_scale=fontsize)
            ForHeatmapPlot = selectMethExp.iloc[:,1:].copy()
            ForHeatmapPlot.columns = ["gene_ID","%s_meth"%(treatment[0]), "%s_exp"%(treatment[0]),"%s_meth"%(treatment[1]),"%s_exp"%(treatment[1]),"Δmethylation","expression","var"]
            p=sns.clustermap(ForHeatmapPlot[["%s_meth"%(treatment[0]),"%s_meth"%(treatment[1]),"%s_exp"%(treatment[0]), "%s_exp"%(treatment[1])]],vmin=0,vmax=100,col_cluster=False,yticklabels=False)
            ForHeatmapPlot.reset_index(inplace=True)
            roworder=p.dendrogram_row.reordered_ind # the order from dendrogram
            redata=selectMethExp.reindex(roworder)
            heatmapplotname = 'heatmap.png'
            plt.savefig(i + "_" + j +"_"+ heatmapplotname, dpi = 300) # save heatmap
            #plt.show()
            print("complete %s on %s. There are %d differential genes"%(context,target,len(selectMethExp)))
            heatmaptablename = "heatmap_rankgenes.txt" 
            redata.iloc[:,1:-1].to_csv(i + "_" + j +"_"+ heatmaptablename, sep="\t",index=False, header=True) # save the order of genes
    return heatmapplotname,heatmaptablename,redata




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
    # get the values of all genes and differential genes
    Samplelist,dropres,selectMeth,selectMethExp, selectMethExp_1, selectMethExp_2, selectMethExp_3, selectMethExp_4, deltaMethmean, deltaExpmean, variation = compare_value(Samplelist=args.samplelist,context=args.context, target=args.target,add=1.0, prob_cutoff = args.prob_cutoff, usecutoff=args.cutoff, Methcutoff=args.meththreshold,Expcutoff= args.expthreshold)
    # create the final heatmap figure
    heatmapplotname, heatmaptablename, selectMethExp = heatmap(Samplelist, dropres, selectMeth, selectMethExp, context=args.context,target=args.target,maxexp=args.emax,maxmeth=args.mmax, Expcutoff=args.expthreshold,fontsize=args.fontsize)
    script_end = time.time() # end time
    print "preprocess time %s"%HowManyTime(script_start,script_end) # total time for analyses




if __name__ == '__main__':
    main()




