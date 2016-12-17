---
layout: page
title: User Guide
subtitle: Everything there is to know
---

# Outline

Basic Principles

* [What does dDocent do?](#what-does-ddocent-do)
	* [Quality Filtering](#quality-filtering)


# What does dDocent do?

## Quality Filtering
dDocent first checks that files are named properly and then uses TrimGalore! to remove low quality bases from reads and remove reads that are primarily adapter sequence.  To perfom quality filtering, simple answer “yes” when dDocent asks “Do you want to quality trim your reads?”  After quality filtering has been performed once, it does not need to be performed again.  Files that have a .R1.fq and .R2.fq are the filtered FASTQ files.





















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
