---
layout: page
title: Download Current Release and Manual Install
subtitle: 
---

# Current Release

Click here to download: <a class="btn btn-success" href="https://github.com/jpuritz/dDocent/archive/v2.2.15.tar.gz"><span class="glyphicon glyphicon-download-alt" aria-hidden="true"></span> dDocent-2.2.15</a>

Alternatively:

```
curl -L -O https://github.com/jpuritz/dDocent/archive/v2.2.15.tar.gz
tar xvzf v2.2.15.tar.gz
```
 
# Requirements
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
|bedtools| http://bedtools.readthedocs.io/en/latest/ |
|vcflib| https://github.com/ekg/vcflib |
|gnuplot| http://www.gnuplot.info |
|gnu-parallel| http://www.gnu.org/software/parallel/ |
|bamtools|https://github.com/pezmaster31/bamtools|
|java| http://www.oracle.com/technetwork/java/javase/downloads/index.html|
|PEAR read merger| http://sco.h-its.org/exelixis/web/software/pear/ |


Also, FreeBayes requires cmake for compiling.  Make sure it is installed on your system. http://www.cmake.org/cmake/resources/software.html
**PEAR neads to be installed as pearRM in your $PATH

# Manual installation.

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
