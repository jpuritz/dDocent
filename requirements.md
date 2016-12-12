---
layout: page
title: Requirements
subtitle: THESE HAVE CHANGED AS OF VERSION 2.0
---

Instead of reinventing the wheel, dDocent relies almost entirely on third party software to complete every step of the 
analysis pipeline, and users are encouraged to familiarize themselves with several of these programs, especially Rainbow, 
BWA, FreeBayes, GATK, and VCFtools.  Below is a list of all the dependencies of dDocent and websites to reference the software:

| Software        | Link                             |
| ------------- |------------------------------------|
|FreeBayes      | https://github.com/ekg/freebayes   |
|STACKS         | http://creskolab.uoregon.edu/stacks|
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
