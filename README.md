# MethGET

A bioinformatics pipeline for correlating DNA methylation and gene expression


# Contents
[System Requirements](#SystemRequirements)
[Installation](#Installation)
[Single-methylome analyses]()
  - Correlation analyses of genome-wide DNA methylation and gene expression
  - Ordinal association analyses with genes ranked by gene expression level
  - Distribution of DNA methylation by groups of genes with different expression level
  - Average methylation level profiles with different expression groups around genes
 

Multiple-methylome analyses:
  - Gene level association between changes of DNA methylation and changes of gene expression
  - Heatmap representation of DNA methylation and gene expression data together
 

# <a name="SystemRequirements"></a>System Requirements
* Linux environment
* CPU：No special restrictions, but CPU has 16 cores is more efficient

* MEM：12GB or higer (for plant sample) / 256GB or higher (for human sample)
* GCC 5.4.0 +

* Python 2.7 
(Type " python -V" to see the installed version. Python2 could be downloaded from  http://www.python.org/download/) 
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

# MethGET Input
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

# Running MethGET
## Data Preprocessing
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
## Single-methylome analyses
#### Correlation analyses of genome-wide DNA methylation and gene expression
#### scatter.py
**Usage:**
```
```
**Example:**
```
```
#### Ordinal association analyses with genes ranked by gene expression level
#### rankscatter.py
**Usage:**
```
```
**Example:**
```
```
#### Distribution of DNA methylation by groups of genes with different expression level
#### grouping.py
**Usage:**
```
```
**Example:**
```
```
#### Average methylation level profiles with different expression groups around genes
#### metaplot.py
**Usage:**
```
```
**Example:**
```
```
## Multiple-methylome analyses
#### Gene level association between changes of DNA methylation and changes of gene expression
#### comparison.py -p scatter
**Usage:**
```
```
**Example:**
```
```
#### Heatmap representation of DNA methylation and gene expression data together
#### comparison.py -p heatmap
**Usage:**
```
```
**Example:**
```
```



