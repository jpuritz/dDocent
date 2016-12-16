---
layout: page
title: User Guide
subtitle: Everything there is to know
---

#Running
If dDocent is installed to your $PATH, change to the data directory and type:

	dDocent

Otherwise it can be run like any other BASH script:

	bash /PATH_TO_dDOCENT/dDocent
#Running with configuration file
The file can be named anything, but must follow the format below:
```bash
Number of Processors
24
Maximum Memory
0
Trimming
no
Assembly?
yes
Type_of_Assembly
PE
Clustering_Similarity%
0.86
Mapping_Reads?
yes
Mapping_Match_Value
1
Mapping_MisMatch_Value
3
Mapping_GapOpen_Penalty
5
Calling_SNPs?
yes
Email
jpuritz@gmail.com
```
Run:
```bash
dDocent config.file
```
#User Guide

For a detailed user guide please see: http://ddocent.wordpress.com
