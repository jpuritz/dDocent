---
layout: page
title: Easy local Install
subtitle: Using bioconda
---


[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/ddocent/README.html)

### Easy local install with bioconda

Install Miniconda: http://conda.pydata.org/miniconda.html

Add the bioconda channel:

```
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda
```

Create a dDocent conda environment:

```
conda create -n ddocent_env ddocent=2.2.4
```

Activate the dDocent environment:

```
source activate ddocent_env
```

Run dDocent:

```
dDocent
```

Close the environment when you're done:

```
source deactivate
```
