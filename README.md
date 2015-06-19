#dDocent version 2.0 has arrived with major updates:

1.  The pipeline now employs a two-step cutoff for data to be included in assembly.
2.  Assembly accuracy has been improved by replacing the sparse seed clustering of rainbow with alignment based clustering in CD-hit
3.  dDocent can now natively handle single end data and paired end data with substantial overlap between paired reads.
4.  Parallelization of variant calling has been improved and is now faster with smaller memory loads.
5.  dDocent can now be run non-interactively by loading a configuration file.

dDocent
=======

This script serves as an interactive bash wrapper to QC, assemble, map, and call SNPs from double digest RAD data.  It is designed to run on Linux based machines with large memory capacity and multiple processing cores

##Check out the tutorials!
https://github.com/jpuritz/dDocent/tree/master/tutorials

#Requirements
**THESE HAVE CHANGED AS OF VERSION 2.0**

Instead of reinventing the wheel, dDocent relies almost entirely on third party software to complete every step of the 
analysis pipeline, and users are encouraged to familiarize themselves with several of these programs, especially Rainbow, 
BWA, FreeBayes, GATK, and VCFtools.  Below is a list of all the dependencies of dDocent and websites to reference the software:

| Software        | Link                             |
| ------------- |------------------------------------|
|FreeBayes      | https://github.com/ekg/freebayes   |
|STACKS         | http://creskolab.uoregon.edu/stacks|
|FastQC		      | http://www.bioinformatics.babraham.ac.uk/projects/fastqc/ |
|Trimmomatic	  | http://www.usadellab.org/cms/?page=trimmomatic |
|Mawk			      | http://invisible-island.net/mawk/ |
|BWA		  	    | http://bio-bwa.sourceforge.net |
|SAMtools		    | http://samtools.sourceforge.net |
|VCFtools		    | http://vcftools.sourceforge.net/index.html |
|rainbow		    | http://sourceforge.net/projects/bio-rainbow/files/ |
|seqtk			    | https://github.com/lh3/seqtk |
|CD-HIT		      | http://weizhong-lab.ucsd.edu/cd-hit/ |
|bedtools| https://code.google.com/p/bedtools/ |
|vcflib| https://github.com/ekg/vcflib |
|gnuplot| http://www.gnuplot.info |
|gnu-parallel| http://www.gnu.org/software/parallel/ |
|bamtools|https://github.com/pezmaster31/bamtools|
|java| http://www.oracle.com/technetwork/java/javase/downloads/index.html|
|PEAR read merger**| http://sco.h-its.org/exelixis/web/software/pear/ |

Also, FreeBayes requires cmake for compiling.  Make sure it is installed on your system. http://www.cmake.org/cmake/resources/software.html
**PEAR neads to be installed as pearRM in your $PATH

#Installation

dDocent is designed to run on a multicore, high memory capacity linux based computer.  As stated above, dDocent depends on several other software packages and assumes that they will be installed in your $PATH directory and that all dDocent dependencies are in a single directory.  The easiest way to do this, for all users of your machine is to install everything into the /usr/local/bin directory.  You will need administrator or "root" privileges to do this.

If you don't have access to the /usr/local/bin directory, don't worry.  dDocent can be installed locally in your user account.  To do this, follow these simple commands:

	cd ~

	mkdir dDocent

	nano .bash_profile

If this file is blank, type:

	PATH="~/dDocent:${PATH}"
	export PATH

Otherwise, simply add ~/dDocent to the end of the existing string.

**Now if you are using a Mac computer, things get a little trickier.  You need to make sure you have Xcode installed, as well as the command line tools.  After this is complete, download the gcc complier from (http://hpc.sourceforge.net) and install it according to the website's instructions.   You will also have to install git from (http://git-scm.com/download).**

If you want more information on setting your $PATH and this setup process, check out the Palumbi Lab’s Simple Fool’s Guide for a good explanation and tutorial on what $PATH is and how to set it (http://sfg.stanford.edu/computer.html).

Once $PATH is setup, there is a VERY simplistic installation script located in the GitHub Repository called install_dDocent_requirements.  To run it, simply type:

	bash install_dDocent_requirements <your path directory>

The script will check to see if any of the required packages are installed and if they aren’t download and install them.  If you are installing computer wide, you probably will need to run the script as sudo.

If all went well, typing “dDocent” and hitting return should start the pipeline.

dDocent requires that your raw data are split up by tagged individual and follow the naming convenction of:

	Pop1_Sample1.F.fq.gz Pop1_Sample1.R.fq.gz

dDocent uses raw reads for reference assembly and trimmed reads for read mapping and SNP/variant calling.  If the user is not using dDocent for trimming, trimmed reads must already be in the directory and must follow the naming convention below:

	Pop1_001.R1.fq.gz  Pop1_001.R2.fq.gz

	Pop1_002.R1.fq.gz  Pop1_002.R2.fq.gz

Where R1 are trimmed forward reads and R2 are trimmed paired-end reads.


These files must all be in the same directory.

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
