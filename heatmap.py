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
    group1.add_argument("-t","--target",type=str,default="Promoter",choices=["Promoter", "Gene_Body","Exon","Intron"],help="choose the target region of methylation, default is Promoter")
    group2 = parser.add_argument_group('Important general arguments')
    group2.add_argument("-mthr","--meththreshold",default="auto",help="set cutoff of differential methylated genes. default 'auto' uses Δ methylation, CG:10, CHG:1, CHH:1")
    group2.add_argument("-ethr","--expthreshold",default=1.0,type=float,help="set cutoff to identify genes that have expression changes, default uses Δ expression log2FC is 1")
    group2.add_argument("-mmax","--mmax",default=100,type=float,help="set the max methylation value for heatmap, default is 100")
    group2.add_argument("-emax","--emax",default=20,type=float,help="set the max expression value for heatmap, default is 20")
    #group2.add_argument("-ad","--addGEvalue",default=1.0,type=float,help="add a small value on gene expression value to calculate the log2FC, default is 1")
    group4 = parser.add_argument_group('Graphing arguments')
    group4.add_argument("--fontsize",default=1.2,type=int,help="fontsize")
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



def heatmap(Samplelist, dropres, selectMeth, selectMethExp, context="CG",target="Gene_Body",maxexp=20 ,maxmeth="auto", Expcutoff=2,fontsize=1.2):
    con=[]
    con.append(context)
    tar=[]
    tar.append(target)
    value={}
    for i in con:
        for j in tar:
            treatment = Samplelist[3].unique()
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
                ## to truncate
                selectMethExp["%sexpmean"%(treatment[0])].loc[selectMethExp["%sexpmean"%(treatment[0])]>maxexp]=maxexp
                selectMethExp["%sexpmean"%(treatment[1])].loc[selectMethExp["%sexpmean"%(treatment[1])]>maxexp]=maxexp
                ###  to let the exp value from (0,20) to (0,100), right! multiply the scale
                selectMethExp["%sexpmean"%(treatment[0])]=selectMethExp["%sexpmean"%(treatment[0])]*(100.0/maxexp)
                selectMethExp["%sexpmean"%(treatment[1])]=selectMethExp["%sexpmean"%(treatment[1])]*(100.0/maxexp)
                ########## methylation
                ## to truncate
                selectMethExp["%smethmean"%(treatment[0])].loc[selectMethExp["%smethmean"%(treatment[0])]>maxmeth]=maxmeth
                selectMethExp["%smethmean"%(treatment[1])].loc[selectMethExp["%smethmean"%(treatment[1])]>maxmeth]=maxmeth
                ###  to let the meth value from (0,20) to (0,100), right! multiply the scale
                selectMethExp["%smethmean"%(treatment[0])]=selectMethExp["%smethmean"%(treatment[0])]*(100.0/maxmeth)
                selectMethExp["%smethmean"%(treatment[1])]=selectMethExp["%smethmean"%(treatment[1])]*(100.0/maxmeth)
            sns.set(font_scale=fontsize)
            ForHeatmapPlot = selectMethExp.iloc[:,1:].copy()
            ForHeatmapPlot.columns = ["gene_ID","ctrl_meth", "ctrl_exp","trt_meth","trt_exp","Δmethylation","expression","var"]
            p=sns.clustermap(ForHeatmapPlot[["ctrl_meth","trt_meth","ctrl_exp", "trt_exp" ]],vmin=0,vmax=100,col_cluster=False,yticklabels=False)
            ForHeatmapPlot.reset_index(inplace=True)
            roworder=p.dendrogram_row.reordered_ind
            redata=selectMeth.reindex(roworder)
            heatmapplotname = 'heatmap.png'
            plt.savefig(i + "_" + j +"_"+ heatmapplotname, dpi = 300)
            #plt.show()
            print("complete %s on %s. There are %d differential genes"%(context,target,len(selectMethExp)))
            heatmaptablename = "heatmap_rankgenes.txt" 
            redata.iloc[:,1:-1].to_csv(i + "_" + j +"_"+ heatmaptablename, sep="\t",index=False, header=True)
            return heatmapplotname,heatmaptablename,redata

def main():
    parser = get_parser()
    args = parser.parse_args()
    Samplelist,dropres,selectMeth,selectMethExp, selectMethExp_1, selectMethExp_2, selectMethExp_3, selectMethExp_4, deltaMethmean, deltaExpmean, variation = compare_value(Samplelist=args.samplelist,context=args.context, target=args.target, add=1.0, Methcutoff=args.meththreshold,Expcutoff= args.expthreshold)
    heatmapplotname, heatmaptablename, selectMethExp = heatmap(Samplelist, dropres, selectMeth, selectMethExp, context=args.context,target=args.target,maxexp=args.emax,maxmeth=args.mmax, Expcutoff=args.expthreshold,fontsize=args.fontsize)





if __name__ == '__main__':
    main()



