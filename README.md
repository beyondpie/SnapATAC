
# SnapATAC Fix Version

## Dependencies
* Use python 2.7 version snaptools

- Though we can install snaptools on python >= 3.0, snaptools have a binary string habit when encoding/decoding the chromosmes, which lead later mapping could be wrong if we ignore this.

- When use `addBmatToSnap`, python2-version snaptools and python3-version snaptools give different colname names. But the matrix is the same. (Only the column order is different, so the matrix is then completely wrong even we align the colnames).

## Install on Linux
Use conda to maintain the enviroment.

```
## SnapATAC depends on raster, which depends on terra, and terra needs gdal.
conda install -c conda-forge libgdal
## install SnapATAC
Rscript -e 'remotes::install_github(repo = "beyondpie/SnapATAC", ref = "szu")'
```

## Install on MacOS
Personnally, I use Brew to maintain the enviroment. I found conda
version R may have some issues on MacOS, so later I just skip it. 

```
brew install gdal proj sqlite
## when install terra, I also have "libproj or sqlite3 not found in
## standard or given location" errors.
## From https://github.com/r-spatial/sf/issues/1312
## In R
install.packages(“terra”, type = "source", configure.args = "--with-gdal-config=/usr/local/bin/gdal-config --with-geos-config=/usr/local/bin/geos-config --with-proj-include=/usr/local/include/ --with-proj-lib=/usr/local/lib/ --with-proj-share=/usr/local/share/proj/")
## Maybe this also works: install.packages(“terra”, configure.args = "--with-proj-lib=/usr/local/lib/")

## Then we can install SnapATAC in R 
remotes::install_github(repo = "beyondpie/SnapATAC", ref = "szu")
```

## Original README:
## SnapATAC (Latest Updates: 2019-09-19)
**SnapATAC** (**S**ingle **N**ucleus **A**nalysis **P**ipeline for **ATAC**-seq) is a fast, accurate and comprehensive method for analyzing single cell ATAC-seq datasets. 

## Latest News
* [SnapATAC links distal elements to putative target genes](https://github.com/r3fang/SnapATAC/blob/master/examples/10X_PBMC_15K/README.md#gene_peak_pair)
* [SnapATAC integrates scRNA and scATAC](https://github.com/r3fang/SnapATAC/blob/master/examples/10X_PBMC_15K/README.md)
* [SnapATAC employs a new method for dimensionality reduction](https://github.com/r3fang/SnapATAC/blob/master/examples/10X_brain_5k/README.md#diffusion_maps)
* [SnapATAC enables clustering using leiden algorithm](https://github.com/r3fang/SnapATAC/blob/master/examples/10X_brain_5k/README.md#cluster)
* [SnapATAC enables batch effect correction](https://github.com/r3fang/SnapATAC/blob/master/examples/10X_snATAC/README.md)
* [SnapATAC enables motif analysis using chromVAR](https://github.com/r3fang/SnapATAC/blob/master/examples/10X_brain_5k/README.md#homer_chromVAR)

## FAQs
* [How to run SnapATAC on 10X dataset?](https://github.com/r3fang/SnapATAC/wiki/FAQs#10X_snap)
* [I already ran CellRanger, can I use its output for SnapATAC?](https://github.com/r3fang/SnapATAC/wiki/FAQs#cellranger_output)
* [How can I analyze combine multiple samples together?](https://github.com/r3fang/SnapATAC/wiki/FAQs#multi_snap)
* [How to group reads from any subset of cells?](https://github.com/r3fang/SnapATAC/wiki/FAQs#group_reads)
* [What is a snap file?](https://github.com/r3fang/SnapATAC/wiki/FAQs#whatissnap)
* [How to create a snap file from fastq file?](https://github.com/r3fang/SnapATAC/wiki/FAQs#CEMBA_snap)

## Requirements  
* Linux/Unix
* Python (>= 2.7 & < 3.0) (SnapTools) (highly recommanded for 2.7);
* R (>= 3.4.0 & < 3.6.0) (SnapATAC) (3.6 does not work for rhdf5 package);

## Pre-print  
Rongxin Fang, Sebastian Preissl, Xiaomeng Hou, Jacinta Lucero, Xinxin Wang, Amir Motamedi, Andrew K. Shiau, Eran A. Mukamel, Yanxiao Zhang, M. Margarita Behrens, Joseph Ecker, Bing Ren. *Fast and Accurate Clustering of Single Cell Epigenomes Reveals Cis-Regulatory Elements in Rare Cell Types.* bioRxiv 615179; doi: https://doi.org/10.1101/615179

## Installation

SnapATAC has two components: [Snaptools](https://github.com/r3fang/SnapTools) and [SnapATAC](https://github.com/r3fang/SnapATAC). 

* SnapTools - a python module for pre-processing and working with [snap](https://github.com/r3fang/SnapATAC/wiki/FAQs) file. 
* SnapATAC  - a R package for the clustering, annotation, motif discovery and downstream analysis.    

Install snaptools from PyPI. See how to install snaptools on [FAQs](https://github.com/r3fang/SnapATAC/wiki/FAQs). 
**NOTE:** Please use python 2.7 if possible. 

```bash
$ pip install snaptools
```

Install SnapATAC R pakcage (development version). 

```
$ R
> library(devtools)
> install_github("r3fang/SnapATAC")
```

## Galleries & Tutorials (click on the image for details)
[<img src="./images/10X_brain_5k.png" width="280" height="318" />](./examples/10X_brain_5k/README.md)
[<img src="./images/PBMC_ATAC_RNA.png" width="280" height="318" />](./examples/10X_PBMC_15K/README.md)
[<img src="./images/10X_snATAC.png" width="280" height="318" />](./examples/10X_snATAC/README.md)
