# MethGET

A bioinformatics pipeline for correlating DNA methylation and gene expression


# Contents
## [System Requirements](#SystemRequirements)

## [Installation](#Installation)

## [MethGET Input](#MethGETInput)

## [Running MethGET](#RunningMethGET)

### [Data Preprocessing](#DataPreprocessing)

### [Single-methylome analyses](#Singlemeth):
  - [Correlation analyses of genome-wide DNA methylation and gene expression](#scatter)
  - [Ordinal association analyses with genes ranked by gene expression level](#rankscatter)
  - [Distribution of DNA methylation by groups of genes with different expression level](#boxplot)
  - [Average methylation level profiles with different expression groups around genes](#metaplot)
 

### [Multiple-methylome analyses](#Multiplemeth):
  - [Gene level association between changes of DNA methylation and changes of gene expression](#deltascatter)
  - [Heatmap representation of DNA methylation and gene expression data together](#heatmap)
 

# <a name="SystemRequirements"></a>System Requirements
* Linux environment

* Python 2.7 
(Type " python -V" to see the installed version. Python2.7 could be downloaded from  http://www.python.org/download/) 
* [CGmaptools](https://cgmaptools.github.io/quick-start/)
* [ucsc-gtftogenepred](https://bioconda.github.io/recipes/ucsc-gtftogenepred/README.html)

Python modules such as "numpy", "pandas", "matplotlib", and "pyBigWig" is needed. To install the packages, use the following commands on linux terminal:
```
	$ pip install numpy
	$ pip install pandas
	$ pip install matplolib
	$ pip install math
	$ pip install scipy
	$ pip install argparse
	$ pip install pyBigWig
```
# <a name="Installation"></a>Installation

# <a name="MethGETInput"></a>MethGET Input
* For single methylome
> 1. DNA Methylation
>>  CGmap.gz file
```
chr1    G       13538   CG      CG      0.6     6       10
chr1    G       13539   CHG     CC      0.0     0       9
chr1    G       13541   CHH     CA      0.0     0       9
chr1    G       13545   CHH     CA      0.0     0       8
```
> 2. Gene Expression File
>> Tab-delimited text ﬁle
```
AT1G01310       0.152868
AT1G01320       9.06088
AT1G01340       14.1157
AT1G01350       24.8099
AT1G01355       0.082233
AT1G01360       42.5391
AT1G01370       8.28111
```
> 3. Gene Annotation File
>> gene annotation in GTF file 

* For multiple methylomes 
>1. Sample Description File 
>> Sample list
```
WT_A    WT_A_CGmap.gz   WT_A_exp.txt    WT
WT_B    WT_B_CGmap.gz   WT_B_exp.txt    WT
MT_A    MT_A_CGmap.gz   MT_A_exp.txt    MT
MT_B    MT_B_CGmap.gz   MT_B_exp.txt    MT
```
>2. Gene Annotation File
>> gene annotation in GTF file


# <a name="RunningMethGET"></a>Running MethGET 
## <a name="DataPreprocessing"></a>Data Preprocessing
#### preprocess.py
Use the script preprocess.py to preprocess the data for downstream analyses

**Usage:**
preprocess.py {-s <samplelist> | -n <samplename> -f <cgmap> -e <expressionfile>} -g <gtf> [options]
```
optional arguments:
  -h, --help            show this help message and exit

For samplelist:
  -s SAMPLELIST, --samplelist SAMPLELIST
                        input sample description file

For individual data:
  -n SAMPLENAME, --samplename SAMPLENAME
                        Setting the name of the set of data. This determines
                        the names of output files for downstream analyses
  -f CGMAPFILENAME, --cgmapfilename CGMAPFILENAME
                        input CGmap file
  -e EXPRESSIONFILE, --expressionfile EXPRESSIONFILE
                        input gene expression file

General options:
  -g GTFFILE, --gtffile GTFFILE
                        input gene annotation file
  -c CUTOFF, --cutoff CUTOFF
                        minimum cytosines that are covered by reads
```
**Example:**
```
# individual data
python preprocess.py -n demo -f WT.CGmap -e WT.exp -g genes.gtf
# for samplelist
python preprocess.py -s samplelist.txt -g genes.gtf
```
## <a name="Singlemeth"></a>Single-methylome analyses
#### <a name="scatter"></a>Correlation analyses of genome-wide DNA methylation and gene expression
#### scatter.py
**Usage:**
```
usage:  [-h] [-n SAMPLENAME] [-p {scatter,kernel}] [-c {CG,CHG,CHH}]
        [-t {Promoter,Gene_Body,Exon,Intron}] [-cor {False,pearson,spearman}]
        [-re0 {True,False}] [-thrs THRESHOLD] [-xlim XLIMIT] [-ylim YLIMIT]
        [--dotsize DOTSIZE] [--textsize TEXTSIZE] [--ticksize TICKSIZE]
        [--labelsize LABELSIZE] [--titlesize TITLESIZE]

optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  -n SAMPLENAME, --samplename SAMPLENAME
                        the name of the set of data
  -p {scatter,kernel}, --plot {scatter,kernel}
                        create scatterplot or kernel density plot, default is
                        'scatterplot'
  -c {CG,CHG,CHH}, --context {CG,CHG,CHH}
                        choose the context of methylation, default is CG
  -t {Promoter,Gene_Body,Exon,Intron}, --target {Promoter,Gene_Body,Exon,Intron}
                        choose the target region of methylation, default is
                        'Promoter'

Important general arguments:
  -cor {False,pearson,spearman}, --correlation {False,pearson,spearman}
                        select the type of correlation, default is pearson
  -re0 {True,False}, --skip0 {True,False}
                        Whether genes with 0 expression value would be
                        included. Default is to include them
  -thrs THRESHOLD, --threshold THRESHOLD
                        Whether skip genes with expression value that is too
                        high, default is to skip genes higher than 1000
                        expression. If want to include them, please set 'None'
  -xlim XLIMIT, --xlimit XLIMIT
                        Nemeric zoom in the gene expression value to clearly
                        understand the distribution
  -ylim YLIMIT, --ylimit YLIMIT
                        Nemeric zoom in the DNA methylation level to clearly
                        understand the distribution

Graphing arguments:
  --dotsize DOTSIZE     dotsize, default is 30
  --textsize TEXTSIZE   textsize, default is 20
  --ticksize TICKSIZE   ticksize, default is 15
  --labelsize LABELSIZE
                        labelsize, default is 20
  --titlesize TITLESIZE
                        titlesize, default is 20
```
**Example:**
```
```
#### <a name="rankscatter"></a>Ordinal association analyses with genes ranked by gene expression level
#### rankscatter.py
**Usage:**
```
usage:  [-h] [-n SAMPLENAME] [-p {scatter,kernel}] [-c {CG,CHG,CHH}]
        [-t {Gene_Body,Promoter,Exon,Intron}] [-re0 {True,False}]
        [-thrs THRESHOLD] [-shsca {True,False}] [-line {True,False}]
        [-smoo SMOOTH_N] [-ylim YLIMIT] [--dotsize DOTSIZE]
        [--textsize TEXTSIZE] [--ticksize TICKSIZE] [--labelsize LABELSIZE]
        [--titlesize TITLESIZE] [--legendsize LEGENDSIZE]

optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  -n SAMPLENAME, --samplename SAMPLENAME
                        the name of the set of data
  -p {scatter,kernel}, --plot {scatter,kernel}
                        create scatterplot or kernel density plot, default is
                        'scatterplot'
  -c {CG,CHG,CHH}, --context {CG,CHG,CHH}
                        choose the context of methylation, default is 'CG'
  -t {Gene_Body,Promoter,Exon,Intron}, --target {Gene_Body,Promoter,Exon,Intron}
                        choose the target region of methylation, default is
                        'Promoter'

Important general arguments:
  -re0 {True,False}, --skip0 {True,False}
                        Whether genes with 0 expression value would be
                        included. Default is to include them
  -thrs THRESHOLD, --threshold THRESHOLD
                        Whether to skip genes with expression value that is
                        too high, default is to skip genes higher than 1000
                        expression. If want to include them, please set 'None'

chart visulaization arguments:
  -shsca {True,False}, --showscatterplot {True,False}
                        whether to show the scatterplot or not, default is to
                        show
  -line {True,False}, --smoothline {True,False}
                        whether to show the fitting curve of methylation and
                        gene expression, default is to show
  -smoo SMOOTH_N, --smooth_N SMOOTH_N
                        set the number of ticks to average when drawing the
                        fitting curve, default is 500
  -ylim YLIMIT, --ylimit YLIMIT
                        Nemeric zoom in the DNA methylation level to clearly
                        understand the distribution

Graphing arguments:
  --dotsize DOTSIZE     dotsize
  --textsize TEXTSIZE   textsize
  --ticksize TICKSIZE   ticksize
  --labelsize LABELSIZE
                        labelsize
  --titlesize TITLESIZE
                        titlesize
  --legendsize LEGENDSIZE
                        legendsize
```
**Example:**
```
```
#### <a name="boxplot"></a>Distribution of DNA methylation by groups of genes with different expression level
#### grouping.py
**Usage:**
```
```
**Example:**
```
```
#### <a name="metaplot"></a>Average methylation level profiles with different expression groups around genes
#### metaplot.py
**Usage:**
```
```
**Example:**
```
```
## <a name="Multiplemeth"></a>Multiple-methylome analyses
#### <a name="deltascatter"></a>Gene level association between changes of DNA methylation and changes of gene expression
#### comparison.py -p scatter
**Usage:**
```
```
**Example:**
```
```
#### <a name="heatmap"></a>Heatmap representation of DNA methylation and gene expression data together
#### comparison.py -p heatmap
**Usage:**
```
```
**Example:**
```
```



