dDocent
=======

This script serves as an interactive bash wrapper to QC, assemble, map, and call SNPs from double digest RAD data.
Instead of reinventing the wheel, dDocent relies almost entirely on third party software to complete every step of the 
analysis pipeline, and users are encouraged to familiarize themselves with several of these programs, especially Rainbow, 
BWA, GATK, and VCFtools.  Below is a list of all the dependencies of dDocent and websites to reference the software:

GATK 			      http://www.broadinstitute.org
STACKS	  	    http://creskolab.uoregon.edu/stacks/
cutadapt		    http://code.google.com/p/cutadapt/
FastQC		      http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
TrimGalore!	    http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/
Mawk			      http://invisible-island.net/mawk/
BWA		  	      http://bio-bwa.sourceforge.net
SAMtools		    http://samtools.sourceforge.net
Picard			    http://picard.sourceforge.net
VCFtools		    http://vcftools.sourceforge.net/index.html
rainbow		      http://sourceforge.net/projects/bio-rainbow/files/
seqtk			      https://github.com/lh3/seqtk
CD-HIT		      http://weizhong-lab.ucsd.edu/cd-hit/
Seq_filter.pl		https://code.google.com/p/seq-filter/downloads/list
cutseq_fasta.pl http://code.google.com/p/nash-bioinformatics-codelets/


dDocent requires that your raw data are split up by tagged individual and follow the naming convenction of:

Pop1_Sample1.F.fq and Pop1_Sample1.R.fq

These files must all be in the same directory.

To run dDocent, simply switch to the directory and type dDocent and press <enter>.


For a detailed user guide please see: http://ddocent.wordpress.com
