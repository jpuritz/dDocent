![logo](logo.png)
_RADSeq Bioinformatics and Beyond_

[![alt text](https://anaconda.org/bioconda/ddocent/badges/downloads.svg)](https://anaconda.org/bioconda/ddocent) 
[![dDocent documentation](https://img.shields.io/badge/documentation-website-informational?logo=Read%20The%20Docs&logoColor=white)](https://www.ddocent.com)

dDocent is simple bash wrapper to QC, assemble, map, and call SNPs from almost any kind of RAD sequencing. If you have a reference already, dDocent can be used to call SNPs from almost any type of NGS data set. It is designed to run on Linux based machines with large memory capacity and multiple processing cores, and it can be modified for use on HPC. 

## Installing

### with `conda` (recommended)
```bash
> conda install -c bioconda dDocent

# or into a fresh environment #

> conda create -n environmentname -c bioconda dDocent

```

### manually
```bash
> chmod +x ./install_dDocent_requirements

> sh ./install_dDocent_requirements
```

## How does dDocent compare?

![ngs comparison](https://github.com/jpuritz/dDocent/blob/master/Sample%20Comparsion.png)

-----

_The "d" is silent_ 🤫
