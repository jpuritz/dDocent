---
layout: page
title: Bioconda Install
subtitle: less than 10 lines of code
---

### The beauty of Conda

Conda is an open source package and environment management system for installing multiple versions of software packages and their dependencies and switching easily between them. It works on Linux, OS X and Windows, and was created for Python programs but can package and distribute any software.

Miniconda is a small version that includes only conda, Python, and the packages they depend on. Over **720** scientific packages and their dependencies can be installed individually from the Continuum repository with the “conda install” command.


### So easy, your grandfather could do this...

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
conda create -n ddocent_env ddocent=2.2.7
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
