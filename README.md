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

Programs with * are only required for dDocent.GATK

#####Dependencies can be installed with the install_dDocent.FB_requirements script by calling it and passing it an installation directory.  Remember to run as sudo if installing to a system directory (for all users).  Also, FreeBayes requires cmake for compiling.  Make sure it is installed on your system. http://www.cmake.org/cmake/resources/software.html


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
