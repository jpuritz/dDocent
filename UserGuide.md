---
layout: page
title: User Guide
subtitle: Everything there is to know
---

# Outline

* [Basic Principles: What does dDocent do?](#what-does-ddocent-do)
	* [Quality Filtering](#quality-filtering)
	* [Assembly](#de-novo-assembly)


# What does dDocent do?

## Quality Filtering
dDocent first checks that files are named properly and then uses TrimGalore! to remove low quality bases from reads and remove reads that are primarily adapter sequence.  To perfom quality filtering, simple answer “yes” when dDocent asks “Do you want to quality trim your reads?”  After quality filtering has been performed once, it does not need to be performed again.  Files that have a .R1.fq and .R2.fq are the filtered FASTQ files.

## *De Novo* Assembly
dDocent uses a novel hybrid approach to assemble reference contigs (RAD loci) based on some of the unique feature of ddRAD fragments.  All reads are concatenated into single forward and reverse FASTA files.  From here, forward and reverse reads are tabulated together, and the complete set of unique reads (loci) are tabulated, along with the number of occurrences of each unique locus.  Reads that do not have a high number of occurrences are likely to be either sequence errors or polymorphisms that are shared by only a few individuals.  This distribution usually follows the asymptotic relationship seen in the figure below:

Screen Shot 2013-10-19 at 1.45.57 AMWhere there are a large proportion of reads that only have one or two occurences, meaning they will not likely be informative on the population scale.  Even if a locus was polymorphic to the point that it was unique to the individual, one can expect that one unique allele would be present at least the level of coverage expected from the sequencing strategy, for example 10X or 20X coverage.  This is where we expect to see the slope of the distribution even out.  In the above example, 10X would be a good coverage to start with, but experimentation for the best value will be needed for each taxa.  Set the cutoff too low, and extraneous reads will be included further down the pipeline and eat up valuable computational time and make it more likely that sequencing errors are included in the data.  Set the cutoff too high, and usuable polymorphic loci may be excluded from subsequent analyses.  NOTE: This data is most appropriate for RAD methods that result in fixed length loci, such as ddRAD or ezRAD.


















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
