dDocent
=======

This script serves as an interactive bash wrapper to QC, assemble, map, and call SNPs from double digest RAD data.

#Requirements

Instead of reinventing the wheel, dDocent relies almost entirely on third party software to complete every step of the 
analysis pipeline, and users are encouraged to familiarize themselves with several of these programs, especially Rainbow, 
BWA, GATK, and VCFtools.  Below is a list of all the dependencies of dDocent and websites to reference the software:

| Software        | Link                             |
| ------------- |------------------------------------|
|GATK           | http://www.broadinstitute.org      |
|STACKS         | http://creskolab.uoregon.edu/stacks|
|cutadapt       | http://code.google.com/p/cutadapt/ |
|FastQC		      | http://www.bioinformatics.babraham.ac.uk/projects/fastqc/ |
|TrimGalore!	  | http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/ |
|Mawk			      | http://invisible-island.net/mawk/ |
|BWA		  	    | http://bio-bwa.sourceforge.net |
|SAMtools		    | http://samtools.sourceforge.net |
|Picard			    | http://picard.sourceforge.net |
|VCFtools		    | http://vcftools.sourceforge.net/index.html |
|rainbow		    | http://sourceforge.net/projects/bio-rainbow/files/ |
|seqtk			    | https://github.com/lh3/seqtk |
|CD-HIT		      | http://weizhong-lab.ucsd.edu/cd-hit/ |
|Seq_filter.pl  | https://code.google.com/p/seq-filter/downloads/list |
|cutseq_fasta.pl| http://code.google.com/p/nash-bioinformatics-codelets/ |


#####Dependencies can be installed with the install_dDocent_requirements script by calling it and passing it an installation directory


dDocent requires that your raw data are split up by tagged individual and follow the naming convenction of:

	Pop1_Sample1.F.fq Pop1_Sample1.R.fq

These files must all be in the same directory.

#Running
If dDocent is installed to your $PATH, change to the data directory and type:

	ddocent 

Otherwise it can be run like any other BASH script:

	sh ./PATH_TO_dDOCENT/dDocent

#User Guide

For a detailed user guide please see: http://ddocent.wordpress.com
