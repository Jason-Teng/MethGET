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
warnings.filterwarnings('ignore')



def get_parser():
    """
    Create a parser and add arguments
    """
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group('Required arguments')
    group1.add_argument("-n","--samplename",type=str,help="the name of the set of data")
    group1.add_argument("-p","--plot",type=str,default="region",choices=["region","site"],help="choose the format of the metaplot, default is region")
    group1.add_argument("-ma","--metavaluefile",default="False",choices=["True","False"],help="put in the metaplot value file if you have created, default is False")  
    group1.add_argument("-nb","--numberofgroup",default=5,type=int,help="define how many group to seperate gene expression, default is 5")
    group1.add_argument("-re0","--skip0",default="False",choices=["True","False"],help="remove genes that their expression value is equal to 0, default is not removing them.")    
    group2 = parser.add_argument_group('Important arguments')
    group2.add_argument("-No_bins","--numberofbins",type=int,help="for 'region' or 'site' plot. 'region' define the total of windows from upstream to downstream, suggest 30. 'site' is the windows of one side, suggest 10.")
    group2.add_argument("-psn","--posname",default="TSS",choices=["TSS","TES"],help="for 'site' plot. Choose the specific position of the metagene, default is TSS (transcription start site)")
    group2.add_argument("-bp","--basepair",default=2000,type=int,help="for 'site' plot. The basepairs flanking by the chosen position, default is 2000")
    group2.add_argument("-yaxis","--yaxissetting",default="auto",choices=["auto","set100"],help="choose if the ticks of yaxis set on 100, default 'auto' will automatically adjusted.")    
    group3 = parser.add_argument_group('Graphing arguments')
    group3.add_argument("-xtick","--xticksize",default=20,type=int,help="the size of the xticks (upstream, gene, downstream)")
    group3.add_argument("-ytick","--yticksize",default=15,type=int,help="the size of the yticks (methylation level)")
    group3.add_argument("-label","--labelsize",default=20,type=int,help="labelsize")
    group3.add_argument("-title","--titlesize",default=20,type=int,help="titlesize")
    group3.add_argument("-legend","--legendsize",default=15,type=int,help="legendsize")
    return parser



def bw_all(name):#(bw_CG,bw_CHG,bw_CHH):
    bw_CG=pyBigWig.open("%s_CG.bw"%(name))
    bw_CHG=pyBigWig.open("%s_CHG.bw"%(name))
    bw_CHH=pyBigWig.open("%s_CHH.bw"%(name))
    bwall={}
    bwall["CG"]=bw_CG  
    bwall["CHG"]=bw_CHG
    bwall["CHH"]=bw_CHH
    return bwall,name




def readAllFinalValue(name):
    matrixinfo=pd.read_csv("%smetavaluefile.txt"%(name),sep="\t",header=None)
    n=len(matrixinfo[0].unique())-1
    All_final_value=np.zeros(len(matrixinfo[0].unique()))
    All_final_value=list(All_final_value)
    for gr in range(len(matrixinfo[0].unique())):
        All_final_value[gr]={}
        for con in ["CG","CHG","CHH"]:
            info=matrixinfo.loc[(matrixinfo[0]==matrixinfo[0].unique()[gr]) & (matrixinfo[1]==con)]
            info=info.values.tolist()[0][2:]
            All_final_value[gr][con]=info
    return All_final_value




def AverageMethGeneBW(bwall,chr,left,right):
    """
    get target methylation site from CGmap
    """
    Meth_value_mean={}
    for i in ["CG","CHG","CHH"]:
        bw=bwall["%s"%(i)]
        Meth_mean=np.nanmean(bw.values(chr,left-1,right)) ##left -1 can calculate containing the left side
        Meth_value_mean["%s"%(i)]=Meth_mean
    return Meth_value_mean 



def Metagenebw(bwall,name,chr,sta,end,dir,No_bins):
    """
    calculate average methylation of each bin
    """
    a=No_bins
    if dir=="+":
        length=end-sta 
        upstream=sta-length*0.5
        downstream=end+length*0.5 
        bins=float((downstream - upstream)/a)#bin's length
        All_position=np.zeros(a+1) 
        All_position[0]=upstream 
        for k in range(1,a+1):
            All_position[k]=All_position[k-1]+bins
        left=All_position[0]
        right=All_position[a]
    else:
        length=end-sta #gene length
        upstream=end+length*0.5
        downstream=sta-length*0.5
        bins=float((upstream - downstream)/a)#bin's length
        All_position=np.zeros(a+1) 
        All_position[0]=upstream 
        for k in range(1,a+1):
            All_position[k]=All_position[k-1]-bins
        left=All_position[a]
        right=All_position[0]
    Meth_value_mean={}
    for con in ["CG","CHG","CHH"]:
        Meth_value_mean["%s"%(con)]=np.zeros(a)
        if dir=="+":
            for i in range(0,a):  
                Meth_value_mean["%s"%(con)][i]=AverageMethGeneBW(bwall,chr,int(All_position[i])-1,int(All_position[i+1])-1)["%s"%(con)]  
        else:
            for i in range(0,a):
                Meth_value_mean["%s"%(con)][i]=AverageMethGeneBW(bwall,chr,int(All_position[i+1]),int(All_position[i]))["%s"%(con)] 
    return Meth_value_mean



def All_value_metageneplot(bwall,name,read_bed,n=5,No_bins=30,skip0=False):    
    No_bins = No_bins*2
    order=["1st","2nd","3rd"]
    for i in range(4,n+1):
        order.append("%dth"%(i)) 
    groupname=order[:n]   
    no_exp=read_bed.loc[read_bed["RPKM"]==0] 
    other_exp=read_bed.loc[read_bed["RPKM"]!=0] 
    sort_other_exp=other_exp.sort_values(by='RPKM')
    sepvalue=[0]  
    groups={} 
    for k in range(n):
        sepvalue.append(sort_other_exp.iloc[(len(sort_other_exp)*(k+1)/n)-1,4]) 
        groups[order[k]]=sort_other_exp.loc[(sort_other_exp["RPKM"]>sepvalue[k]) & (sort_other_exp["RPKM"]<=sepvalue[k+1])]
    """
    create final value
    """
    All_final_value=[]
    if skip0==False:
        groupname.insert(0,"no")
        no_exp_All_MethLevel={}
        no_exp_Final_value={} 
        for con in ["CG","CHG","CHH"]:
            no_exp_All_MethLevel["%s"%(con)]=np.zeros((len(no_exp),No_bins))
            for j in range(1,len(no_exp)+1):    
                try:
                    no_exp_All_MethLevel["%s"%(con)][j-1,:]=Metagenebw(bwall,name,no_exp.iloc[j-1,0],no_exp.iloc[j-1,1],no_exp.iloc[j-1,2],no_exp.iloc[j-1,5],No_bins)["%s"%(con)]
                    #print no_exp.iloc[j-1,:], j
                except:
                    pass
                    #print(no_exp.iloc[j-1,:], "error@no", j)
            no_exp_Final_value["%s"%(con)]=np.zeros(No_bins)
            no_exp_All_MethLevel["%s"%(con)]=np.ma.masked_array(no_exp_All_MethLevel["%s"%(con)],np.isnan(no_exp_All_MethLevel["%s"%(con)]))
            for m in range(0,No_bins):
                no_exp_Final_value["%s"%(con)][m]=no_exp_All_MethLevel["%s"%(con)][:,m].mean()*100
        All_final_value.append(no_exp_Final_value)
    for i in range(n):  
        other_exp_All_MethLevel={}
        other_exp_Final_value={}
        for con in ["CG","CHG","CHH"]:
            other_exp_All_MethLevel["%s"%(con)]=np.zeros((len(groups[order[i]]),No_bins))
            for j in range(1,len(groups[order[i]])+1):  
                try:
                    other_exp_All_MethLevel["%s"%(con)][j-1,:]=Metagenebw(bwall,name,groups[order[i]].iloc[j-1,0],groups[order[i]].iloc[j-1,1],groups[order[i]].iloc[j-1,2],groups[order[i]].iloc[j-1,5],60)["%s"%(con)]
                    #print groups[order[i]].iloc[j-1,:], j
                except:
                    pass
                    #print(groups[order[i]].iloc[j-1,:], "error@others", j)
            other_exp_Final_value["%s"%(con)]=np.zeros(No_bins)
            other_exp_All_MethLevel["%s"%(con)]=np.ma.masked_array(other_exp_All_MethLevel["%s"%(con)],np.isnan(other_exp_All_MethLevel["%s"%(con)]))
            for m in range(0,No_bins):
                other_exp_Final_value["%s"%(con)][m]=other_exp_All_MethLevel["%s"%(con)][:,m].mean()*100
        All_final_value.append(other_exp_Final_value)
    matrixinfo=pd.DataFrame()
    matrix = []
    for gr in range(len(groupname)):
        for con in ["CG","CHG","CHH"]:
            info = list(All_final_value[gr][con])
            info.insert(0,"%s"%(con))
            info.insert(0,"%s"%(groupname[gr]))
            matrix.append(info)
    matrixinfo = matrixinfo.append(matrix)
    matrixinfo.to_csv("%smetavaluefile.txt"%(name),header=False,index=False,sep="\t")
    return All_final_value
    




def metageneplot(All_final_value,name,n=5,No_bins=30,skip0=False,yaxis="auto",figuresize=(8,6),xticksize=20,yticksize=15,labelsize=20,titlesize=20,legendsize=15):
    No_bins = No_bins*2
    order=["1st","2nd","3rd"]
    for i in range(4,n+1):
        order.append("%dth"%(i))  
    groupname=order[:n]
    if skip0==False:
        groupname.insert(0,"no")   
    for con in ["CG","CHG","CHH"]:
        maxall=[]
        for j in range(len(All_final_value)): 
            All_final_value[j]
            maxall.append(max(All_final_value[j]["%s"%(con)]))
        highest_value=max(maxall)
        plt.figure(figsize=figuresize)
        plt.subplot(111)
        ### ylim(set the value of y lim) ###
        ylimit=highest_value*3/2
        if ylimit>100:  
            plt.ylim((0, 100))
        else:
            plt.ylim((0, ylimit))
        l1=plt.axvline(x=No_bins*0.25, color='lightgray',linestyle='--')
        l2=plt.axvline(x=No_bins*0.75, color='lightgray',linestyle='--')
        color=cm.rainbow(np.linspace(0.6,0.0,len(All_final_value)))
        color = color.tolist()
        for j in range(len(All_final_value)):
            plt.plot(All_final_value[j]["%s"%(con)],color=color[j],label=groupname[j])
        plt.legend(loc='upper right',fontsize=legendsize)
        plt.title("%s"%(con),fontsize=titlesize)
        plt.xticks([No_bins*0.125, No_bins*0.5,No_bins*0.875 ],[r'$upstream$', r'$gene$', r'$downstream$'],fontsize=xticksize) 
        plt.yticks(fontsize=yticksize)
        plt.ylabel('Methylation Level (%)',fontsize=labelsize)
        plt.savefig('%s_meta_gene_plot_%s_%dexpress.png'%(name,con,n),dpi=300)
        #plt.show()




def Metapointbw(bwall,chr,pos,dir,bp,No_bins):
    totalBin=2*No_bins 
    if dir=="+":
        upstream=pos-bp#upstream position
        downstream=pos+bp #downstream position
        bins=float(bp/No_bins)#bin's length
        All_position=np.zeros(totalBin+1) 
        All_position[0]=upstream
        for k in range(1,totalBin+1):
            All_position[k]=All_position[k-1]+bins
        left=All_position[0] 
        right=All_position[totalBin]
    else:
        #length=end-sta #gene length
        upstream=pos+bp
        downstream=pos-bp
        bins=float(bp/No_bins) #bin's length
        All_position=np.zeros(totalBin+1) 
        All_position[0]=downstream
        for k in range(1,totalBin+1): 
            All_position[k]=All_position[k-1]-bins
        left=All_position[totalBin]
        right=All_position[0]
    Meth_value_mean={} 
    for con in ["CG","CHG","CHH"]: 
        Meth_value_mean["%s"%(con)]=np.zeros(totalBin)
        if dir=="+":
            for i in range(0,totalBin):
                Meth_value_mean["%s"%(con)][i]=AverageMethGeneBW(bwall,chr,int(All_position[i])-1,int(All_position[i+1])-1)["%s"%(con)]
        else:
            for i in range(0,totalBin):
                Meth_value_mean["%s"%(con)][i]=AverageMethGeneBW(bwall,chr,int(All_position[i+1]),int(All_position[i]))["%s"%(con)]
    return Meth_value_mean






def All_value_metapointplot(bwall,name,read_bed,n=5,posname="TSS",bp=2000,No_bins=10,skip0=False):
    totalBin=2*No_bins
    order=["1st","2nd","3rd"]
    for i in range(4,n+1):
        order.append("%dth"%(i))  
    groupname=order[:n]
    no_exp=read_bed.loc[read_bed["RPKM"]==0]
    other_exp=read_bed.loc[read_bed["RPKM"]!=0]
    sort_other_exp=other_exp.sort_values(by='RPKM')
    sepvalue=[0]
    groups={}
    for k in range(n):
        sepvalue.append(sort_other_exp.iloc[(len(sort_other_exp)*(k+1)/n)-1,4])
        groups[order[k]]=sort_other_exp.loc[(sort_other_exp["RPKM"]>sepvalue[k]) & (sort_other_exp["RPKM"]<=sepvalue[k+1])]
    """
    create final value
    """
    All_final_value=[]
    if skip0==False:
        groupname.insert(0,"no")
        no_exp_All_MethLevel={}
        no_exp_Final_value={}
        for con in ["CG","CHG","CHH"]:
            no_exp_All_MethLevel["%s"%(con)]=np.zeros((len(no_exp),totalBin))
            for j in range(1,len(no_exp)+1):    
                if posname=="TSS":
                    pos=no_exp.iloc[j-1,1]
                if posname=="TES":
                    pos=no_exp.iloc[j-1,2]
                try:
                    no_exp_All_MethLevel["%s"%(con)][j-1,:]=Metapointbw(bwall,no_exp.iloc[j-1,0],pos,no_exp.iloc[j-1,5],bp,No_bins)["%s"%(con)]
                    #print no_exp.iloc[j-1,:], j
                except:
                    pass
                    #print no_exp.iloc[j-1,:], "error@no", j
            no_exp_Final_value["%s"%(con)]=np.zeros(totalBin)
            #print no_exp_Final_value["%s"%(con)]
            no_exp_All_MethLevel["%s"%(con)]=np.ma.masked_array(no_exp_All_MethLevel["%s"%(con)],np.isnan(no_exp_All_MethLevel["%s"%(con)]))
            for m in range(0,totalBin):
                no_exp_Final_value["%s"%(con)][m]=no_exp_All_MethLevel["%s"%(con)][:,m].mean()*100
        All_final_value.append(no_exp_Final_value)
    for i in range(n):  
        other_exp_All_MethLevel={}
        other_exp_Final_value={}
        for con in ["CG","CHG","CHH"]:
            other_exp_All_MethLevel["%s"%(con)]=np.zeros((len(groups[order[i]]),totalBin))
            for j in range(1,len(groups[order[i]])+1):  
                if posname=="TSS":
                    pos=groups[order[i]].iloc[j-1,1]
                if posname=="TES":
                    pos=groups[order[i]].iloc[j-1,2]
                try:
                    other_exp_All_MethLevel["%s"%(con)][j-1,:]=Metapointbw(bwall,groups[order[i]].iloc[j-1,0],pos,groups[order[i]].iloc[j-1,5],bp,No_bins)["%s"%(con)]
                    #print other_exp.iloc[j-1,:], j,i
                except:
                    pass
                    #groups[order[i]].iloc[j-1,:], "error@others", j
            other_exp_Final_value["%s"%(con)]=np.zeros(totalBin)
            other_exp_All_MethLevel["%s"%(con)]=np.ma.masked_array(other_exp_All_MethLevel["%s"%(con)],np.isnan(other_exp_All_MethLevel["%s"%(con)]))
            for m in range(0,totalBin):
                other_exp_Final_value["%s"%(con)][m]=other_exp_All_MethLevel["%s"%(con)][:,m].mean()*100
        All_final_value.append(other_exp_Final_value)
    matrixinfo=pd.DataFrame()
    matrix = []
    for gr in range(len(groupname)):
        for con in ["CG","CHG","CHH"]:
            info = list(All_final_value[gr][con])
            info.insert(0,"%s"%(con))
            info.insert(0,"%s"%(groupname[gr]))
            matrix.append(info)
    matrixinfo = matrixinfo.append(matrix)
    matrixinfo.to_csv("%smetavaluefile.txt"%(name),header=False,index=False,sep="\t")
    return All_final_value




def metapointplot(All_final_value, name, posname="TSS",n=5,bp=2000,No_bins=10,skip0=False,yaxis="auto",figuresize=(8,6),xticksize=20,yticksize=15,labelsize=20,titlesize=20,legendsize=15):
    order=["1st","2nd","3rd"]
    for i in range(4,n+1):
        order.append("%dth"%(i))
    groupname=order[:n]
    if skip0==False:
        groupname.insert(0,"no")
    for i in ["CG","CHG","CHH"]:
        #### max ####
        maxall=[]
        for j in range(len(All_final_value)):  
            All_final_value[j]
            maxall.append(max(All_final_value[j]["%s"%(i)]))
        highest_value=max(maxall)
        plt.figure(figsize=figuresize)
        plt.subplot(111)
        ### ylim(set the value of y lim) ###
        ylimit=highest_value*3/2
        if ylimit>100:
            plt.ylim((0, 100))
        else:
            plt.ylim((0, ylimit))
        l1=plt.axvline(x=No_bins-0.5, color='lightgray',linestyle='--')
        color=cm.rainbow(np.linspace(0.6,0.0,n+1))
        color = color.tolist()
        for j in range(len(All_final_value)):
            plt.plot(All_final_value[j]["%s"%(i)],color=color[j],label=groupname[j])
        plt.legend(loc='upper right',fontsize=legendsize)
        plt.title("%s"%(i),fontsize=titlesize)
        plt.xticks([0,No_bins-0.5,2*No_bins-1],["-%d"%(bp),'%s'%(posname),"+%d"%(bp)],fontsize=xticksize)
        plt.yticks(fontsize=yticksize)
        plt.ylabel('Methylation Level (%)',fontsize=labelsize)
        plt.savefig('%s_meta_%s_plot_%s_%dexpress.png'%(name,posname,i,n),dpi = 300)
        #plt.show()





def HowManyTime(tbegin,tend):
    tTotal=tend-tbegin
    tsec=tTotal%60
    ttolmin=tTotal//60
    thour=ttolmin//60
    tmin=ttolmin%60
    suretime="running time is %d hour, %d minutes, %.2f seconds"%(thour,tmin,tsec)
    return suretime



def main():
    tbegin = time.time() 
    parser = get_parser()
    args = parser.parse_args()
    exp=pd.read_csv("%s_exp.txt"%(args.samplename),sep="\t",names=["gene_ID","RPKM"])
    read_bed=pd.read_csv('gb.bed', sep='\t',names=["chr","left","right","gene_ID","RPKM","direction"],dtype={"chr":str})
    res=pd.merge(read_bed[["chr","left","right","gene_ID","direction"]],exp,on=["gene_ID"]) ##how='inner'
    res=res[["chr","left","right","gene_ID","RPKM","direction"]]
    bwall,name=bw_all(args.samplename)
    if args.plot == "region":
        if args.metavaluefile=="False":
            All_final_value=All_value_metageneplot(bwall,name,read_bed=res,n=args.numberofgroup,No_bins=args.numberofbins,skip0=False)
            if eval(args.skip0) ==True:
                no_group = All_final_value.pop(0)
        if args.metavaluefile=="True":
            All_final_value=readAllFinalValue(name)
            if eval(args.skip0) ==True:
                no_group = All_final_value.pop(0)
        metageneplot(All_final_value,name,n=args.numberofgroup,No_bins=args.numberofbins,skip0=eval(args.skip0),yaxis="auto",figuresize=(8,6),xticksize=args.xticksize,yticksize=args.yticksize,labelsize=args.labelsize,titlesize=args.labelsize,legendsize=args.legendsize)
    if args.plot == "site":
        if args.metavaluefile=="False":
            All_final_value=All_value_metapointplot(bwall,name,read_bed=res,n=args.numberofgroup, posname=args.posname, bp=args.basepair, No_bins=args.numberofbins,skip0=False)
            if eval(args.skip0) ==True:
                no_group = All_final_value.pop(0)
        if args.metavaluefile=="True":
            All_final_value=readAllFinalValue(name)
            if eval(args.skip0) ==True:
                no_group = All_final_value.pop(0)
        metapointplot(All_final_value,name,n=args.numberofgroup,posname=args.posname, bp=args.basepair,No_bins=args.numberofbins,skip0=eval(args.skip0),yaxis="auto",figuresize=(8,6),xticksize=args.xticksize,yticksize=args.yticksize,labelsize=args.labelsize,titlesize=args.labelsize,legendsize=args.legendsize)
    tend = time.time()
    #print(HowManyTime(tbegin,tend))



if __name__ == '__main__':
    main()


