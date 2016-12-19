---
layout: page
title: User Guide
subtitle: Everything there is to know
bigimg: /ddlogo.png
---

# Outline

* [Basic Principles: What does dDocent do?](#what-does-ddocent-do)
	* [Quality Filtering](#quality-filtering)
	* [Assembly](#de-novo-assembly)
	* [Mapping](#read-mapping)
	* [SNP calling](#snp-calling)
	* [SNP Filtering](#snp-filtering)
* [Getting Started](#getting-started)
	* [Requirements](#requirements)
		* [Raw Sequences](#raw-sequences)
		* [Naming convention](#naming-convention)
	* [Running](#running-dDocent)
* [Customization](#Customizing-dDocent)

# What does dDocent do?

## Quality Filtering
dDocent first checks that files are named properly and then uses TrimGalore! to remove low quality bases from reads and remove reads that are primarily adapter sequence.  To perfom quality filtering, simple answer “yes” when dDocent asks “Do you want to quality trim your reads?”  After quality filtering has been performed once, it does not need to be performed again.  Files that have a .R1.fq and .R2.fq are the filtered FASTQ files.

## *De Novo* Assembly
dDocent uses a novel hybrid approach to assemble reference contigs (RAD loci) based on some of the unique feature of ddRAD fragments.  All reads are concatenated into single forward and reverse FASTA files.  From here, forward and reverse reads are tabulated together, and the complete set of unique reads (loci) are tabulated, along with the number of occurrences of each unique locus.  Reads that do not have a high number of occurrences are likely to be either sequence errors or polymorphisms that are shared by only a few individuals.  This distribution usually follows the asymptotic relationship seen in the figure below:

Screen Shot 2013-10-19 at 1.45.57 AMWhere there are a large proportion of reads that only have one or two occurences, meaning they will not likely be informative on the population scale.  Even if a locus was polymorphic to the point that it was unique to the individual, one can expect that one unique allele would be present at least the level of coverage expected from the sequencing strategy, for example 10X or 20X coverage.  This is where we expect to see the slope of the distribution even out.  In the above example, 10X would be a good coverage to start with, but experimentation for the best value will be needed for each taxa.  Set the cutoff too low, and extraneous reads will be included further down the pipeline and eat up valuable computational time and make it more likely that sequencing errors are included in the data.  Set the cutoff too high, and usuable polymorphic loci may be excluded from subsequent analyses.  NOTE: This data is most appropriate for RAD methods that result in fixed length loci, such as ddRAD or ezRAD.

After the user sets the data cutoff, dDocent removes the extraneous sequence reads, splits the reads back into forward and reverse pairs and then inputs them into the program rainbow using its default parameters.  Rainbow cluster forward reads based on similarity (6 bp difference with default setting).  These clusters are then recursively divided based on reverse reads into groups representing single alleles.  Reads in merged contigs are then assembled using a greedy algorithm.  The longest contig for each cluster is finally selected as the representative reference sequence for that contig.  If the forward read does not overlap with the reverse read (almost always the case with ddRAD), the forward read is pasted to the reverse read with a ten N basepairs as padding.  Finally, reference sequences are clustered based on overall sequence similarity using the program CD-HIT.  Alternatively, de novo assembly can be skipped and the user can provide a fasta file with reference sequences.

## Read Mapping

Reads are mapped (aligned to reference sequences) using the MEM algorithm of BWA.  The user can choose to set alternative values for the match score value, mismatch score, and gap opening penalty.  According the the BWA manual, the sequence error rate is approximately: {.75 * exp[-log(4) * B/A]}.  Again, experimentation for finding the optimal values for a particular taxa is recommended.  The default settings are meant for mapping reads to the human genome and are fairly conservative.  BWA outputs SAM files which are then converted to BAM files using SAMtools.

## SNP Calling

dDocent uses a scatter gather technique to speed up SNP/INDEL calling.  In short, reference contigs are split up into a single file for each processing core.  FreeBayes then calls variants for each genomic interval simultaneously (changes to default parameters include setting a minimum mapping and base quality score to PHRED 10). FreeBayes is a Bayesian-based variant detection program that uses assembled haplotype sequences to simultaneously call SNPs, INDELS, multi-nucleotide polymorphisms (MNPs), and complex events (composite insertion and substitution events) from alignment files; FreeBayes has the added benefit for population genomics of using reads across multiple individuals to improve genotyping.  After all instances of FreeBayes finish, raw SNP/INDEL calls are concatenated into a single variant call file (VCF) using VCFtools.

## SNP Filtering

Final SNP data sets will depend on the individual project, so dDocent does only minimal filtering.  Using VCFtools, SNPs are filtered to only those that are called in 90% of all individuals.  They can be found in Final.recode.vcf; however, these are mainly for diagnostic purposes between runs.  VCFtools can be used to filter SNP and INDEL calls using a variety of criteria and it is recommended that users familiarize themselves with the program to produce a truly final call set.

For more ideas on stringent filtering, I recommend checking out the [SNP Filtering Tutorial](/filtering)


# Getting Started

This section of the user guide will walk you through using the pipeline.

## Requirements

After you have dDocent properly installed (see [Bioconda Install](/bioconda) and [Manual Install](/manual)), you need to create a working directory:

```bash
mkdir my_dDocent_working_dir
```

### Raw Sequences
In this directory, you need to place **RAW** and **DEMULTIPLEXED** sequencing files.  

**Trimmed reads will fail *de novo* assembly** 
If performing *de novo* assembly, it's essential that no read trimming or adapter removal has taken place before the dDocent pipeline.  If a reference is being supplied, then trimmed reads may be used.

### Naming Convention
dDocent requires that sequence files be gzipped FASTQ format and the files **MUST MUST MUST** follow a specific naming convention.File names must contain a locality/population identifier and an individual identifier, and these two identifiers must be separated by a single `_`

For example:

```
Pop1_Sample1.F.fq.gz Pop1_Sample1.R.fq.gz
```

Please don't use `_` in either of the identifiers.  

dDocent uses raw reads for reference assembly and trimmed reads for read mapping and SNP/variant calling.  If dDocent is not being used for trimming, trimmed reads must already be in the directory and must follow the naming convention below:

	Pop1_001.R1.fq.gz  Pop1_001.R2.fq.gz

	Pop1_002.R1.fq.gz  Pop1_002.R2.fq.gz

Where R1 are trimmed forward reads and R2 are trimmed paired-end reads (If PE sequencing is being used).



# Running

Simply type:

````
dDocent
```
Now, dDocent will start interactive configuration by asking you questions:

```
dDocent 2.2.7 

Contact jpuritz@gmail.com with any problems 

 
Checking for required software

All required software is installed!

dDocent run started Sat Dec 17 23:28:22 EST 2016 

40 individuals are detected. Is this correct? Enter yes or no and press [ENTER]
```

To start, dDocent is checking for all the required software and then counting the number of individuals that it finds in the current directory. This is just a basic check to make sure you're in the correct directory.  Simple enter `yes` to continue

```
dDocent detects 72 processors available on this system.
Please enter the maximum number of processors to use for this analysis.
```

Next, dDocent tries to automatically detect the number of processors available on your system.  It then lets you set a hard limit for the number of processors to use simultaneously during a dDocent run.  Simply enter the number you want to continue.

```
32
```
```
dDocent detects 251G maximum memory available on this system.
Please enter the maximum memory to use for this analysis. The size can be postfixed with 
K, M, G, T, P, k, m, g, t, or p which would multiply the size with 1024, 1048576, 1073741824, 
1099511627776, 1125899906842624, 1000, 1000000, 1000000000, 1000000000000, or 1000000000000000 respectively.
For example, to limit dDocent to ten gigabytes, enter 10G or 10g
```

dDocent tries to automatically detect the maximum amount of memory on your system.  It then lets you set a hard limit for the amount of memory that it will use.  This will only apply to SNP calling.

#### This option does not work on all systems!!!!!
`GNU-Parallel` implements this option and it unfortunately does not work properly on all `LINUX` distributions, and any value other than `0` may cause SNP calling to hold indefinately.  Unless this has ben tested on your system, I recommend entering:

```
0
```







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





Outputs:

dDocent will output several different files as part of the pipeline.  The main outputs of interest are:

TotalRawSNPs.vcf- This file, in the standard Variant Call Format, has the raw SNP, INDEL, MNP, and complex variant calls for every individual.  This is the file that will be used for further filtering (see VCFtools or vcflib) to produce the final data set.  It is important to note that FreeBayes combines SNP and INDEL calls that are in within a default 3bp window into haplotype calls of the complex variant calls.  To properly look at SNPs only, complex variants need to be decomposed with vcfallelicprimatives from the vcflib package (https://github.com/ekg/vcflib) and then INDELs can be filtered with VCFtools or vcflib.

Final.recode.vcf- This is the filtered VCF file from the very end of the pipeline.  It's useful for comparing different runs of the pipeline.

X.bam- These files are Binary Alignment Maps (BAMs).  They are the mapping of reads for every individual to the reference contigs.  BAM files are input directly into FreeBayes for variant calling.  Note, this file cannot be viewed directly; use SAMtools for inspection.

reference.fasta- These are the contigs generated by the de novo assembly in FASTA format.

uniqseq.data- This file is a tab delimited file the level of coverage in the first column and the number of unique sequences with more than that level of coverage is the second column.

For information on the other outputs, refer to the code comments in the dDocent.FB file.  They are  useful for debugging purposes.

Customizing dDocent:

Further customization of the pipeline can be achieved by modifying the script directly.  Any text editor can be used.  Below is a guide on how to modify specific sections of code:

Read Trimming:

By default, dDocent looks for Illumina TruSeq adapters and trims off basepairs with a quality score less than 10.  Both BWA and FreeBayes take base quality scores into account, so excessive trimming is not necessary nor recommended.  To modify this find line 351.  Settings can be modified with the following parameters:

-a          first read adapter sequence

-a2        second read adapter sequence

-q         quality score for trimming

Visit the trim galore website (http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) for more information on settings.

Assembly:

By default, dDocent only alters one parameter from the defaults of rainbow for assembly, and that is to change the maximum number of mismatches to cluster reads together for assembly from 4 to 6.  This setting appears to work well with polymorphic marine species.  Users can change this value on line 506 by changing the number after the -m parameter.  Other parameters can be fine tuned in rainbow if necessary on lines 507 and 508.  Please see the rainbow documentation and paper for more information about fine tuning the assembly.

Read mapping:

Through the command line interface, dDocent users already are able to customize several read mapping parameters for BWA.  BWA does have additional options that can be added to line 237  and line 506 (Changes need to be made to both), though they shouldn't be necessary for most circumstances.  Please consult the BWA manual for more information.

SNP/INDEL calling:

FreeBayes is a highly customizable variant calling program and can be adapted to many different needs.  To see all the different possible options of FreeBayes, simply type:

freebayes -h

This brings up the FreeBayes help menu that lists and explains all options.  Options that may be of interest are listed below:

--populations FILE     If you are expecting your data to be highly structured, you may consider using this option.  FreeBayes does do some population specific calculations for heterozygosity levels and mutation rates.  This allows you to partition the population model.  To do this, you must supply a file that lists the sample and it's population.  Note, you should only use this is you have a strong a priori reason to expect substantial population structure.

-E      This parameter affects how haplotypes are called in the VCF file.  The default value is 3, meaning that variants within 3 bp of each other will be called as a single contiguous haplotype.   This represents the highest tradeoff between haplotype length and sensitivity.  See this discussion for more information.

-m     This parameter controls the minimum mapping quality score for reads to be considered for genotyping.  By default, dDocent sets this to PHRED 10 (or 90% probability of being true).

-q      This sets the minimum base quality score for a bp to be used for genotyping.  By default, dDocent sets this to PHRED 10 (or 90% probability of being true).

-V     This parameter tells FreeBayes to ignore read placement bias and strand bias in the calculation of SNP site quality scores.  By default, dDocent leaves this parameter in because it, at worst, conservatively biases quality scores.  However, this may be too conservative because one would expect extreme strand bias in RAD sequencing.

To customize FreeBayes these parameters can be changed or altered on lines 287  and 297 (These lines of code must be identical).

VCF Filtering:

To customize the basic VCF filtering performed by dDocent, simply edit line 317. Please see the VCFtools documentation for a list of options.
