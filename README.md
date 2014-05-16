dDocent
=======

This script serves as an interactive bash wrapper to QC, assemble, map, and call SNPs from double digest RAD data.  It is designed to run on Linux based machines with large memory capacity and multiple processing cores

There are now two different versions of dDocent: dDocent.FB and dDocent.GATK.  dDocent.FB uses minimal BAM file preparation steps before calling SNPs and INDELS simultaneously using FreeBayes (Garrison & Marth 2012).  dDocent.GATK uses GATK (McKenna et al. 2010) for INDEL realignment, SNP and INDEL genotyping (using HaplotypeCaller), and variant quality score recalibration, largely following GATK Best Practices recommendations (DePristo et al. 2011; Auwera & Carneiro 2013).  The modules represent two different strategies for SNP/INDEL calling, and are completely independent of one another.

For now, I will be focusing on dDocent.FB because it is substantially faster and has less dependecies.  See http://bcbio.wordpress.com/2013/10/21/updated-comparison-of-variant-detection-methods-ensemble-freebayes-and-minimal-bam-preparation-pipelines/ for a great comparison of FreeBayes and GATK.

#Requirements

Instead of reinventing the wheel, dDocent relies almost entirely on third party software to complete every step of the 
analysis pipeline, and users are encouraged to familiarize themselves with several of these programs, especially Rainbow, 
BWA, FreeBayes, GATK, and VCFtools.  Below is a list of all the dependencies of dDocent and websites to reference the software:

| Software        | Link                             |
| ------------- |------------------------------------|
|FreeBayes      | https://github.com/ekg/freebayes   |
|GATK*          | http://www.broadinstitute.org      |
|STACKS         | http://creskolab.uoregon.edu/stacks|
|cutadapt       | http://code.google.com/p/cutadapt/ |
|FastQC		      | http://www.bioinformatics.babraham.ac.uk/projects/fastqc/ |
|TrimGalore!	  | http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/ |
|Mawk			      | http://invisible-island.net/mawk/ |
|BWA		  	    | http://bio-bwa.sourceforge.net |
|SAMtools		    | http://samtools.sourceforge.net |
|Picard*		    | http://picard.sourceforge.net |
|VCFtools		    | http://vcftools.sourceforge.net/index.html |
|rainbow		    | http://sourceforge.net/projects/bio-rainbow/files/ |
|seqtk			    | https://github.com/lh3/seqtk |
|CD-HIT		      | http://weizhong-lab.ucsd.edu/cd-hit/ |
|Seq_filter.pl  | https://code.google.com/p/seq-filter/downloads/list |
|cutseq_fasta.pl| http://code.google.com/p/nash-bioinformatics-codelets/ |
|bedtools| https://code.google.com/p/bedtools/ |
|vcflib| https://github.com/ekg/vcflib |

Programs with * are only required for dDocent.GATK

Also, FreeBayes requires cmake for compiling.  Make sure it is installed on your system. http://www.cmake.org/cmake/resources/software.html

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

	sh install_dDocent.FB_requirements <your path directory>

The script will check to see if any of the required packages are installed and if they aren’t download and install them.  If you are installing computer wide, you probably will need to run the script as sudo.

If all went well, typing “dDocent.FB” and hitting return should start the pipeline.


dDocent requires that your raw data are split up by tagged individual and follow the naming convenction of:

	Pop1_Sample1.F.fq Pop1_Sample1.R.fq

These files must all be in the same directory.

#Running
If dDocent is installed to your $PATH, change to the data directory and type:

	dDocent.FB 

Otherwise it can be run like any other BASH script:

	sh /PATH_TO_dDOCENT/dDocent.FB

#User Guide

For a detailed user guide please see: http://ddocent.wordpress.com
