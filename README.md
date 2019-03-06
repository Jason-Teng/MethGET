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
> 1. DNA Methylatioin
>>  CGmap.gz file
> 2. Gene Expression File
>> tab-delimited text ﬁle
> 3. Gene Annotation File
>> GTF file 

* For multiple methylomes 
>1. Sample Description File 
>> sample list
>2. Gene Annotation File
>> GTF file

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
```
**Example:**
```
```
#### <a name="rankscatter"></a>Ordinal association analyses with genes ranked by gene expression level
#### rankscatter.py
**Usage:**
```
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



