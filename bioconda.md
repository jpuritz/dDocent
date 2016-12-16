---
layout: page
title: Bioconda Install
subtitle: less than 10 lines of code
---

### Procedure

Install Miniconda: [http://conda.pydata.org/miniconda.html](http://conda.pydata.org/miniconda.html)

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
