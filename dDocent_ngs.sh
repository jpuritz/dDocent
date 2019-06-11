#!/usr/bin/env bash
export LC_ALL=en_US.UTF-8

##########dDocent##########
VERSION='2.2.22'
#This script serves as an interactive bash wrapper to QC, assemble, map, and call SNPs from double digest RAD (SE or PE), ezRAD (SE or PE) data, or SE RAD data.
#It requires that your raw data are split up by tagged individual and follow the naming convention of:

#Pop_Sample1.F.fq and Pop_Sample1.R.fq

#Prints out title and contact info
echo -e "dDocent" $VERSION "\n"
echo -e "Contact jpuritz@gmail.com with any problems \n\n "

###Code to check for the required software for dDocent

echo "Checking for required software"
DEP=(freebayes mawk bwa samtools vcftools rainbow gnuplot gawk seqtk cd-hit-est bamToBed bedtools parallel vcfcombine bamtools pearRM)
NUMDEP=0
for i in "${DEP[@]}"
do
	if which $i &> /dev/null; then
		foo=0
	else
    		echo "The dependency" $i "is not installed or is not in your" '$PATH'"."
    		NUMDEP=$((NUMDEP + 1))
	fi
done

if find ${PATH//:/ } -maxdepth 1 -name trimmomatic*jar 2> /dev/null| grep -q 'trim' ; then
	TRIMMOMATIC=$(find ${PATH//:/ } -maxdepth 1 -name trimmomatic*jar 2> /dev/null | head -1)
	else
    echo "The dependency trimmomatic is not installed or is not in your" '$PATH'"."
    NUMDEP=$((NUMDEP + 1))
	fi

if find ${PATH//:/ } -maxdepth 1 -name TruSeq2-PE.fa 2> /dev/null | grep -q 'Tru' ; then
	ADAPTERS=$(find ${PATH//:/ } -maxdepth 1 -name TruSeq2-PE.fa 2> /dev/null | head -1)
	else
    echo "The file listing adapters (included with trimmomatic) is not installed or is not in your" '$PATH'"."
    NUMDEP=$((NUMDEP + 1))
    fi

SAMV1=$(samtools 2>&1 >/dev/null | grep Ver | sed -e 's/Version://' | cut -f2 -d " " | sed -e 's/-.*//' | cut -c1)
SAMV2=$(samtools 2>&1 >/dev/null | grep Ver | sed -e 's/Version://' | cut -f2 -d " " | sed -e 's/-.*//' | cut -c3)
	if [ "$SAMV1"  -ge "1" ]; then
		if [ "$SAMV2"  -lt "3" ]; then
        	echo "The version of Samtools installed in your" '$PATH' "is not optimized for dDocent."
        	echo "Please install at least version 1.3.0"
			echo -en "\007"
			echo -en "\007"
			exit 1
		fi
	
	else
		    echo "The version of Samtools installed in your" '$PATH' "is not optimized for dDocent."
        	echo "Please install at least version 1.3.0"
			echo -en "\007"
			echo -en "\007"
			exit 1
	fi

RAINV=(`rainbow | head -1 | cut -f2 -d' ' `)	
	if [[ "$RAINV" != "2.0.2" && "$RAINV" != "2.0.3" && "$RAINV" != "2.0.4" ]]; then
        	echo "The version of Rainbow installed in your" '$PATH' "is not optimized for dDocent."
        	echo -en "\007"
			echo -en "\007"
			echo -en "\007"
        	echo "Is the version of rainbow installed newer than 2.0.2?  Enter yes or no."
			read TEST
			if [ "$TEST" != "yes" ]; then 
        		echo "Please install a version newer than 2.0.2"
        		exit 1
        	fi
        fi
FREEB=(`freebayes | grep -oh 'v[0-9].*' | cut -f1 -d "." | sed 's/v//' `)	
	if [ "$FREEB" != "1" ]; then
        	echo "The version of FreeBayes installed in your" '$PATH' "is not optimized for dDocent."
        	echo "Please install at least version 1.0.0"
        	exit 1
        fi         	
VCFTV=$(vcftools | grep VCF | grep -oh '[0-9]*[a-z]*)$' | sed 's/[a-z)]//')
	if [ "$VCFTV" -lt "10" ]; then
        	echo "The version of VCFtools installed in your" '$PATH' "is not optimized for dDocent."
        	echo "Please install at least version 0.1.11"
        	exit 1
        elif [ "$VCFTV" == "11" ]; then
                VCFGTFLAG="--geno" 
        elif [ "$VCFTV" -ge "12" ]; then
                VCFGTFLAG="--max-missing"
	fi
BWAV=$(bwa 2>&1 | mawk '/Versi/' | sed 's/Version: //g' | sed 's/0.7.//g' | sed 's/-.*//g' | cut -c 1-2)
	if [ "$BWAV" -lt "13" ]; then
        	echo "The version of bwa installed in your" '$PATH' "is not optimized for dDocent."
        	echo "Please install at least version 0.7.13"
        	exit 1
	fi

BTC=$( bedtools --version | mawk '{print $2}' | sed 's/v//g' | cut -f1,2 -d"." | sed 's/2\.//g' )
	if [ "$BTC" -ge "26" ]; then
		BEDTOOLSFLAG="NEW"
		elif [ "$BTC" == "23" ]; then
		BEDTOOLSFLAG="OLD"
		elif [ "$BTC" != "23" ]; then
		echo "The version of bedtools installed in your" '$PATH' "is not optimized for dDocent."
		echo "Please install version 2.23.0 or version 2.26.0 and above"
		exit 1	
	fi
		
if ! awk --version | fgrep -v GNU &>/dev/null; then
         awk=gawk
    else
         awk=awk
fi


if [ $NUMDEP -gt 0 ]; then
	echo -e "\nPlease install all required software before running dDocent again."
	exit 1
else
	echo -e "\nAll required software is installed!"
fi

#This code checks for individual fastq files follow the correct naming convention and are gziped
TEST=$(ls *.fq 2> /dev/null | wc -l )
if [ "$TEST" -gt 0 ]; then
echo -e "\ndDocent is now configured to work on compressed sequence files.  Please run gzip to compress your files."
echo "This is as simple as 'gzip *.fq'"
echo "Please rerun dDocent after compressing files."
exit 1
fi

#Count number of individuals in current directory
NumInd=$(ls *.F.fq.gz | wc -l)
NumInd=$(($NumInd - 0))

#Create list of sample names
ls *.F.fq.gz > namelist
sed -i'' -e 's/.F.fq.gz//g' namelist
#Create an array of sample names
NUMNAMES=$(mawk '/_/' namelist | wc -l)

if [ "$NUMNAMES" -eq "$NumInd" ]; then
	NAMES=( `cat "namelist" `)
else
	echo "Individuals do not follow the dDocent naming convention."
	echo "Please rename individuals to: Locality_Individual.F.fq.gz"
	echo "For example: LocA_001.F.fq.gz"
	exit 1
fi

#Wrapper for main program functions.  This allows the entire file to be read first before execution
main(){
##########User Input Section##########
#This code gets input from the user and assigns variables
######################################

#Sets a start time variable
STARTTIME=$(date)


echo -e "\ndDocent run started" $STARTTIME "\n"


#dDocent can now accept a configuration file instead of running interactively
#Checks if a configuration file is being used, if not asks for user input
if [ -n "$1" ]; then
	CONFIG=$1
	NUMProc=$(grep -A1 Processor $CONFIG | tail -1)
    	MAXMemory=$(grep -A1 Memory $CONFIG | tail -1)
	TRIM=$(grep -A1 Trim $CONFIG | tail -1)
	ASSEMBLY=$(grep -A1 '^Assembly' $CONFIG | tail -1)
	ATYPE=$(grep -A1 Type $CONFIG | tail -1)
	simC=$(grep -A1 Simi $CONFIG | tail -1)
	MAP=$(grep -A1 Mapping_R $CONFIG | tail -1)
	optA=$(grep -A1 _Match $CONFIG | tail -1)
	optB=$(grep -A1 MisMatch $CONFIG | tail -1)
	optO=$(grep -A1 Gap $CONFIG | tail -1)
	SNP=$(grep -A1 SNP $CONFIG | tail -1)
	MAIL=$(grep -A1 Email $CONFIG | tail -1)
	if [ "$ASSEMBLY" == "no" ]; then
		#Prints instructions on how to move analysis to background and disown process
		echo "At this point, all configuration information has been entered and dDocent may take several hours to run." 
		echo "It is recommended that you move this script to a background operation and disable terminal input and output."
		echo "All data and logfiles will still be recorded."
		echo "To do this:"
		echo "Press control and Z simultaneously"
		echo "Type 'bg' without the quotes and press enter"
		echo "Type 'disown -h' again without the quotes and press enter"
		echo ""
		echo "Now sit back, relax, and wait for your analysis to finish."
	else
		echo "dDocent will require input during the assembly stage.  Please wait until prompt says it is safe to move program to the background."
	fi

else
	GetInfo 
fi

#Creates (or appends to) a dDcoent run file recording variables
echo "Variables used in dDocent Run at" $STARTTIME >> dDocent.runs
echo "Number of Processors" >> dDocent.runs
echo $NUMProc >> dDocent.runs
echo "Maximum Memory" >> dDocent.runs
echo $MAXMemory >> dDocent.runs
echo "Trimming" >> dDocent.runs
echo $TRIM >> dDocent.runs
echo "Assembly?" >> dDocent.runs
echo $ASSEMBLY >> dDocent.runs
echo "Type_of_Assembly" >> dDocent.runs
echo $ATYPE >> dDocent.runs
echo "Clustering_Similarity%" >> dDocent.runs
echo $simC >> dDocent.runs
echo "Mapping_Reads?" >> dDocent.runs
echo $MAP >> dDocent.runs
echo "Mapping_Match_Value" >> dDocent.runs
echo $optA >> dDocent.runs
echo "Mapping_MisMatch_Value" >> dDocent.runs
echo $optB >> dDocent.runs
echo "Mapping_GapOpen_Penalty" >> dDocent.runs
echo $optO >> dDocent.runs
echo "Calling_SNPs?" >> dDocent.runs
echo $SNP >> dDocent.runs
echo "Email" >> dDocent.runs
echo $MAIL >> dDocent.runs


##Section of logic statements that dictates the order and function of processing the pipeline

if [[ "$TRIM" == "yes" && "$ASSEMBLY" == "yes" ]]; then
        echo "Trimming reads and simultaneously assembling reference sequences"        
        TrimReads & 2> trim.log
        Assemble
        #setupRainbow 2> rainbow.log
        wait
fi

if [[ "$TRIM" == "yes" && "$ASSEMBLY" != "yes" ]]; then
        echo "Trimming reads"
        TrimReads 2> trim.log
fi                
                
if [[ "$TRIM" != "yes" && "$ASSEMBLY" == "yes" ]]; then                
        Assemble
        #setupRainbow 2> rainbow.log
fi

#Checks to see if reads will be mapped.
if [ "$MAP" != "no" ]; then
echo "Using BWA to map reads."
	if [ reference.fasta -nt reference.fasta.fai ]; then
        samtools faidx reference.fasta
        bwa index reference.fasta &> index.log
	fi
#dDocent now checks for trimmed read files before attempting mapping
        if [[ "$MAP" != "no" && ! -f "${NAMES[@]:(-1)}".R1.fq.gz ]]; then
        	echo "dDocent cannot locate trimmed reads files"
        	echo "Please rerun dDocent with quality trimming"
        	exit 1
        fi
#This next section of code checks to see if the reference was assembled by dDocent 
#and if so, modifies the expected insert length distribution for BWA's metric for proper pairing
        if head -1 reference.fasta | grep -e 'dDocent' reference.fasta 1>/dev/null; then
        	rm lengths.txt &> /dev/null
        	for i in "${NAMES[@]}";
        		do
        		if [ -f "$i.R.fq.gz" ]; then
        		zcat $i.R.fq.gz | head -2 | tail -1 >> lengths.txt
        		fi
        		done	
        	if [ -f "lengths.txt" ]; then
        	MaxLen=$(mawk '{ print length() | "sort -rn" }' lengths.txt| head -1)
        	INSERT=$(($MaxLen * 2 ))
        	INSERTH=$(($INSERT + 100 ))
        	INSERTL=$(($INSERT - 100 ))
        	SD=$(($INSERT / 5))
        	fi
#BWA for mapping for all samples.  As of version 2.0 can handle SE or PE reads by checking for PE read files
        	for i in "${NAMES[@]}"
        	do
        	if [ -f "$i.R2.fq.gz" ]; then
        		bwa mem reference.fasta $i.R1.fq.gz $i.R2.fq.gz -t $NUMProc -a -M -T 20 -A $optA -B $optB -O $optO -R "@RG\tID:$i\tSM:$i\tPL:Illumina" 2> bwa.$i.log | samtools view -@$NUMProc -q 1 -SbT reference.fasta - > $i.bam 2>$i.bam.log
        	else
        		bwa mem reference.fasta $i.R1.fq.gz -t $NUMProc -a -M -T 20 -A $optA -B $optB -O $optO -R "@RG\tID:$i\tSM:$i\tPL:Illumina" 2> bwa.$i.log | samtools view -@$NUMProc -q 1 -SbT reference.fasta - > $i.bam 2>$i.bam.log
        	fi
        	samtools sort -@$NUMProc $i.bam -o $i.bam 
		mv $i.bam $i-RG.bam
		samtools index $i-RG.bam
        	done
        else
        	for i in "${NAMES[@]}"
        	do
        	if [ -f "$i.R2.fq.gz" ]; then
        		bwa mem reference.fasta $i.R1.fq.gz $i.R2.fq.gz -t $NUMProc -a -M -T 20 -A $optA -B $optB -O $optO -R "@RG\tID:$i\tSM:$i\tPL:Illumina" 2> bwa.$i.log | samtools view -@$NUMProc -q 1 -SbT reference.fasta - > $i.bam 2>$i.bam.log
        	else
        		bwa mem reference.fasta $i.R1.fq.gz -t $NUMProc -a -M -T 20 -A $optA -B $optB -O $optO -R "@RG\tID:$i\tSM:$i\tPL:Illumina" 2> bwa.$i.log | samtools view -@$NUMProc -q 1 -SbT reference.fasta - > $i.bam 2>$i.bam.log
        	fi
        	samtools sort -@$NUMProc $i.bam -o $i.bam 
		mv $i.bam $i-RG.bam
		samtools index $i-RG.bam
        	done
        fi
fi

##Creating mapping intervals if needed, CreateIntervals function is defined later in script
#If mapping is being performed, intervals are created automatically

if [ "$MAP" != "no" ]; then
echo "Creating alignment intervals"
ls *-RG.bam >bamlist.list
CreateIntervals 
fi

##SNP Calling Section of code

if [ "$SNP" != "no" ]; then
	#Create list of BAM files
	ls *-RG.bam >bamlist.list
	#If mapping is not being performed, but intervals do not exist they are created
	if [[ "$MAP" == "no" && ! -f "cat-RRG.bam" ]]; then
		CreateIntervals 
	fi
	#Check for runs from older versions to ensure the recreation of cat-RRG.bam
	if [[ "$MAP" == "no" && -f "map.bed" ]]; then
		CreateIntervals 
	fi
	#Check to make sure interval files have been created
	if [[ "$MAP" == "no" && ! -f "mapped.bed" ]]; then
		bamToBed -i cat-RRG.bam > map.bed
		bedtools merge -i map.bed > mapped.bed
		rm map.bed
	fi
	#This code estimates the coverage of reference intervals and removes intervals in 0.01% of depth
	#This allows genotyping to be more effecient and eliminates extreme copy number loci from the data
	if [ "cat-RRG.bam" -nt "cov.stats" ]; then
		if [ "$BEDTOOLSFLAG" == "OLD" ]; then
			coverageBed -abam cat-RRG.bam -b mapped.bed -counts > cov.stats
		else
			bedtools coverage -b  cat-RRG.bam -a mapped.bed -counts -sorted > cov.stats
		fi
	fi
		
	if head -1 reference.fasta | grep -e 'dDocent' reference.fasta 1>/dev/null; then
	
		DP=$(mawk '{print $4}' cov.stats | sort -rn | perl -e '$d=.001;@l=<>;print $l[int($d*@l)]')
		CC=$( mawk -v x=$DP '$4 < x' cov.stats | mawk '{len=$3-$2;lc=len*$4;tl=tl+lc} END {OFMT = "%.0f";print tl/"'$NUMProc'"}')
	else
		DP=$(mawk '{print $4}' cov.stats | sort -rn | perl -e '$d=.00005;@l=<>;print $l[int($d*@l)]')
		CC=$( mawk -v x=$DP '$4 < x' cov.stats | mawk '{len=$3-$2;lc=len*$4;tl=tl+lc} END {OFMT = "%.0f";print tl/"'$NUMProc'"}')
	fi
	mawk -v x=$DP '$4 < x' cov.stats |sort -V -k1,1 -k2,2 | mawk -v cutoff=$CC 'BEGIN{i=1} 
	{
	len=$3-$2;lc=len*$4;cov = cov + lc
	if ( cov < cutoff) {x="mapped."i".bed";print $1"\t"$2"\t"$3 > x}
	else {i=i+1; x="mapped."i".bed"; print $1"\t"$2"\t"$3 > x; cov=0}
	}' 
	
	FB2=$(( $NUMProc / 4 ))

	echo "Using FreeBayes to call SNPs"

	#Creates a population file to use for more accurate genotype calling
	
	cut -f1 -d "_" namelist > p
	paste namelist p > popmap
	rm p
	


###New implementation of SNP calling here to save on memory	
	call_genos(){
	samtools view -@$FB2 -b -1 -L mapped.$1.bed -o split.$1.bam cat-RRG.bam
	samtools index split.$1.bam
	freebayes -b split.$1.bam -t mapped.$1.bed -v raw.$1.vcf -f reference.fasta -m 5 -q 5 -E 3 --min-repeat-entropy 1 -V --populations popmap -n 10
	rm split.$1.bam*
	}
	export -f call_genos

	ls mapped.*.bed | sed 's/mapped.//g' | sed 's/.bed//g' | shuf | parallel --env call_genos --memfree $MAXMemory -j $NUMProc --no-notice call_genos {}
####	
	#ls mapped.*.bed | sed 's/mapped.//g' | sed 's/.bed//g' | shuf | parallel --memfree $MAXMemory -j $FB1 --no-notice --delay 1 freebayes -L bamlist.list -t mapped.{}.bed -v raw.{}.vcf -f reference.fasta -m 5 -q 5 -E 3 --min-repeat-entropy 1 -V --populations popmap -n 10
	#ls mapped.*.bed | sed 's/mapped.//g' | sed 's/.bed//g' | shuf | parallel --memfree $MAXMemory -j $FB1 --no-notice "samtools view -b -L mapped.{}.bed | freebayes -c -t mapped.{}.bed -v raw.{}.vcf -f reference.fasta -m 5 -q 5 -E 3 --min-repeat-entropy 1 -V --populations popmap -n 10"


	rm mapped.*.bed 

	mv raw.1.vcf raw.01.vcf
	mv raw.2.vcf raw.02.vcf
	mv raw.3.vcf raw.03.vcf
	mv raw.4.vcf raw.04.vcf
	mv raw.5.vcf raw.05.vcf
	mv raw.6.vcf raw.06.vcf
	mv raw.7.vcf raw.07.vcf
	mv raw.8.vcf raw.08.vcf
	mv raw.9.vcf raw.09.vcf

	vcfcombine raw.*.vcf | sed -e 's/	\.\:/	\.\/\.\:/g' > TotalRawSNPs.vcf

	if [ ! -d "raw.vcf" ]; then
		mkdir raw.vcf
	fi

	mv raw.*.vcf ./raw.vcf

	echo "Using VCFtools to parse TotalRawSNPS.vcf for SNPs that are called in at least 90% of individuals"
	vcftools --vcf TotalRawSNPs.vcf $VCFGTFLAG 0.9 --out Final --recode --non-ref-af 0.001 --max-non-ref-af 0.9999 --mac 1 --minQ 30 --recode-INFO-all &>VCFtools.log
fi

##Checking for possible errors

if [ "$MAP" != "no" ]; then
ERROR1=$(mawk '/developer/' bwa* 2>/dev/null | wc -l 2>/dev/null) 
fi
ERROR2=$(mawk '/error/' *.bam.log 2>/dev/null | wc -l 2>/dev/null)
ERRORS=$(($ERROR1 + $ERROR2))

#Move various log files to own directory
if [ ! -d "logfiles" ]; then
mkdir logfiles
fi
mv *.txt *.log log ./logfiles 2> /dev/null

#Sending a completion email

if [ $ERRORS -gt 0 ]; then
        echo -e "dDocent has finished with errors in" `pwd` "\n\ndDocent started" $STARTTIME "\n\ndDocent finished" `date` "\n\nPlease check log files\n\n" `mawk '/After filtering, kept .* out of a possible/' ./logfiles/VCFtools.log` "\n\ndDocent" $VERSION "\nThe 'd' is silent, hillbilly." | mailx -s "dDocent has finished with ERRORS!" $MAIL
else
        echo -e "dDocent has finished with an analysis in" `pwd` "\n\ndDocent started" $STARTTIME "\n\ndDocent finished" `date` "\n\n" `mawk '/After filtering, kept .* out of a possible/' ./logfiles/VCFtools.log` "\n\ndDocent" $VERSION "\nThe 'd' is silent, hillbilly." | mailx -s "dDocent has finished" $MAIL
fi


}

##Function definitions

#Function for trimming reads using trimmomatic
trim_reads(){
	TRIMMOMATIC=$(find ${PATH//:/ } -maxdepth 1 -name trimmomatic*jar 2> /dev/null | head -1)
    ADAPTERS=$(find ${PATH//:/ } -maxdepth 1 -name TruSeq2-PE.fa 2> /dev/null | head -1)

	if [ -f $1.R.fq.gz ]; then	
		java -jar $TRIMMOMATIC PE -threads 2 -phred33 $1.F.fq.gz $1.R.fq.gz $1.R1.fq.gz $1.unpairedF.fq.gz $1.R2.fq.gz $1.unpairedR.fq.gz ILLUMINACLIP:$ADAPTERS:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:5:10 $TW &> $1.trim.log
	else 
		java -jar $TRIMMOMATIC SE -threads 2 -phred33 $1.F.fq.gz $1.R1.fq.gz ILLUMINACLIP:$ADAPTERS:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:5:10 $TW &> $1.trim.log
	fi 
}
	
	export -f trim_reads

TrimReads () { 
	#STACKS adds a strange _1 or _2 character to the end of processed reads, this looks for checks for errant characters and replaces them.
	#This functionality is now parallelized and will run if only SE sequences are used.

	STACKS=$(cat namelist| parallel -j $NUMProc --no-notice "zcat {}.F.fq.gz | head -1" | mawk '$0 !~ /\/1$/ && $0 !~ /\/1[ ,	]/ && $0 !~ / 1:.*[A-Z]*/' | wc -l )
	FB1=$(( $NUMProc / 2 ))
	if [ $STACKS -gt 0 ]; then
		
		echo "Removing the _1 character and replacing with /1 in the name of every sequence"
		cat namelist | parallel -j $FB1 --no-notice "zcat {}.F.fq.gz | sed -e 's:_1$:/1:g' > {}.F.fq"
		rm -f *.F.fq.gz
		cat namelist | parallel -j $FB1 --no-notice "gzip {}.F.fq"
	fi

	if [ -f "${NAMES[@]:(-1)}".R.fq.gz ]; then
	
		STACKS=$(cat namelist| parallel -j $NUMProc --no-notice "zcat {}.R.fq.gz | head -1" | mawk '$0 !~ /\/2$/ && $0 !~ /\/2[ ,	]/ && $0 !~ / 2:.*[A-Z]*/'| wc -l )

		if [ $STACKS -gt 0 ]; then
			echo "Removing the _2 character and replacing with /2 in the name of every sequence"
			cat namelist | parallel -j $FB1 --no-notice "zcat {}.R.fq.gz | sed -e 's:_2$:/2:g' > {}.R.fq"
			rm -f *.R.fq.gz
			cat namelist | parallel -j $FB1 --no-notice "gzip {}.R.fq"
		fi
	fi

	cat namelist | parallel -j $NUMProc "zcat {}.F.fq.gz | head -2 | tail -1 >> lengths.txt"
	MLen=$(mawk '{ print length() | "sort -rn" }' lengths.txt| head -1)
    MLen=$(($MLen / 2))
	TW="MINLEN:$MLen"
	cat namelist | parallel --env trim_reads -j $FB1 trim_reads {}	
	mkdir unpaired &>/dev/null
	mv *unpaired*.gz ./unpaired &>/dev/null	
}


#Main function for assembly
Assemble()
{
AWK1='BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}'
AWK2='!/>/'
AWK3='!/NNN/'
PERLT='while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}'
SED1='s/^[ 	]*//'
SED2='s/ /	/g'
FRL=$(zcat ${NAMES[0]}.F.fq.gz | mawk '{ print length() | "sort -rn" }' | head -1)

if [ ${NAMES[@]:(-1)}.F.fq.gz -nt ${NAMES[@]:(-1)}.uniq.seqs ];then
	if [[ "$ATYPE" == "PE" || "$ATYPE" == "RPE" ]]; then
	#If PE assembly, creates a concatenated file of every unique for each individual in parallel
		cat namelist | parallel --no-notice -j $NUMProc "zcat {}.F.fq.gz | mawk '$AWK1' | mawk '$AWK2' > {}.forward"
		cat namelist | parallel --no-notice -j $NUMProc "zcat {}.R.fq.gz | mawk '$AWK1' | mawk '$AWK2' > {}.reverse"
		if [ "$ATYPE" = "RPE" ]; then
			cat namelist | parallel --no-notice -j $NUMProc "paste -d '-' {}.forward {}.reverse | mawk '$AWK3'| sed 's/-/NNNNNNNNNN/' | sort | uniq -c -w $FRL| sed -e '$SED1' | sed -e '$SED2' > {}.uniq.seqs"
		else
			cat namelist | parallel --no-notice -j $NUMProc "paste -d '-' {}.forward {}.reverse | mawk '$AWK3'| sed 's/-/NNNNNNNNNN/' | perl -e '$PERLT' > {}.uniq.seqs"
		fi
		rm *.forward
		rm *.reverse
	fi
	if [ "$ATYPE" == "SE" ]; then
	#if SE assembly, creates files of every unique read for each individual in parallel
		cat namelist | parallel --no-notice -j $NUMProc "zcat {}.F.fq.gz | mawk '$AWK1' | mawk '$AWK2' | perl -e '$PERLT' > {}.uniq.seqs"
	fi
	if [ "$ATYPE" == "OL" ]; then
	#If OL assembly, dDocent assumes that the marjority of PE reads will overlap, so the software PEAR is used to merge paired reads into single reads
		for i in "${NAMES[@]}";
        		do
        		zcat $i.R.fq.gz | head -2 | tail -1 >> lengths.txt
        		done	
        	MaxLen=$(mawk '{ print length() | "sort -rn" }' lengths.txt| head -1)
		LENGTH=$(( $MaxLen / 3))
		for i in "${NAMES[@]}"
			do
			pearRM -f $i.F.fq.gz -r $i.R.fq.gz -o $i -j $NUMProc -n $LENGTH 
			done
		cat namelist | parallel --no-notice -j $NUMProc "mawk '$AWK1' {}.assembled.fastq | mawk '$AWK2' | perl -e '$PERLT' > {}.uniq.seqs"
	fi
		
	
fi

#Create a data file with the number of unique sequences and the number of occurrences

if [ -f "uniq.seqs.gz" ]; then
	if [ uniq.seqs.gz -nt uniq.seqs ]; then
	gunzip uniq.seqs.gz 2>/dev/null
	fi
fi

if [ ! -f "uniq.seqs" ]; then
	cat *.uniq.seqs > uniq.seqs
fi
	

for i in {2..20};
do 
echo $i >> pfile
done
cat pfile | parallel -j $NUMProc --no-notice "echo -n {}xxx && mawk -v x={} '\$1 >= x' uniq.seqs | wc -l" | mawk  '{gsub("xxx","\t",$0); print;}'| sort -g > uniqseq.data
rm pfile


#Plot graph of above data
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale
set xrange [2:20] 
unset label
set title "Number of Unique Sequences with More than X Coverage (Counted within individuals)"
set xlabel "Coverage"
set ylabel "Number of Unique Sequences"
plot 'uniqseq.data' with lines notitle
pause -1
EOF


echo -en "\007"
echo -en "\007"
echo -en "\007"
echo -e "Please choose data cutoff.  In essence, you are picking a minimum (within individual) coverage level for a read (allele) to be used in the reference assembly"

read CUTOFF
if [ "$ATYPE" == "RPE" ]; then
	parallel --no-notice -j $NUMProc mawk -v x=$CUTOFF \''$1 >= x'\' ::: *.uniq.seqs | cut -f2 | sort | uniq -c -w $FRL | sed -e 's/^[ 	]*//' | sed -e 's/ /	/g' > uniqCperindv
else
	parallel --no-notice -j $NUMProc mawk -v x=$CUTOFF \''$1 >= x'\' ::: *.uniq.seqs | cut -f2 | perl -e 'while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}' > uniqCperindv
fi
if [ "$NumInd" -gt 10 ]; then
	NUM=$(($NumInd / 2))
else
	NUM=$NumInd
fi

for ((i = 2; i <= $NUM; i++));
do
echo $i >> ufile
done

cat ufile | parallel -j $NUMProc --no-notice "echo -n {}xxx && mawk -v x={} '\$1 >= x' uniqCperindv | wc -l" | mawk  '{gsub("xxx","\t",$0); print;}'| sort -g > uniqseq.peri.data
rm ufile


#Plot graph of above data

gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Number of Unique Sequences present in more than X Individuals"
set xlabel "Number of Individuals"
set ylabel "Number of Unique Sequences"
plot 'uniqseq.peri.data' with lines notitle
pause -1
EOF

echo -en "\007"
echo -en "\007"
echo -en "\007"
echo -e "Please choose data cutoff.  Pick point right before the assymptote. A good starting cutoff might be 10% of the total number of individuals"

read CUTOFF2

#Prints instructions on how to move analysis to background and disown process
echo "At this point, all configuration information has been entered and dDocent may take several hours to run." 
echo "It is recommended that you move this script to a background operation and disable terminal input and output."
echo "All data and logfiles will still be recorded."
echo "To do this:"
echo "Press control and Z simultaneously"
echo "Type 'bg' without the quotes and press enter"
echo "Type 'disown -h' again without the quotes and press enter"
echo ""
echo "Now sit back, relax, and wait for your analysis to finish."

#Now that data cutoffs have been chosen, reduce data set to specified set of unique reads, convert to FASTA format,
#and remove reads with substantial amounts of adapters

mawk -v x=$CUTOFF2 '$1 >= x' uniqCperindv > uniq.k.$CUTOFF.c.$CUTOFF2.seqs
cut -f2 uniq.k.$CUTOFF.c.$CUTOFF2.seqs > totaluniqseq
mawk '{c= c + 1; print ">dDocent_Contig_" c "\n" $1}' totaluniqseq > uniq.full.fasta
LENGTH=$(mawk '!/>/' uniq.full.fasta  | mawk '(NR==1||length<shortest){shortest=length} END {print shortest}')
LENGTH=$(($LENGTH * 3 / 4))
$awk 'BEGIN {RS = ">" ; FS = "\n"} NR > 1 {print "@"$1"\n"$2"\n+""\n"gensub(/./, "I", "g", $2)}' uniq.full.fasta > uniq.fq
java -jar $TRIMMOMATIC SE -threads $NUMProc -phred33 uniq.fq uniq.fq1 ILLUMINACLIP:$ADAPTERS:2:30:10 MINLEN:$LENGTH
mawk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' uniq.fq1 > uniq.fasta
mawk '!/>/' uniq.fasta > totaluniqseq
rm uniq.fq*

#If this is a PE assebmle
if [[ "$ATYPE" == "PE" || "$ATYPE" == "RPE" ]]; then
	#Reads are first clustered using only the Forward reads using CD-hit instead of rainbow
	sed -e 's/NNNNNNNNNN/	/g' uniq.fasta | cut -f1 > uniq.F.fasta
	CDHIT=$(python -c "print(max("$simC" - 0.1,0.8))")
	cd-hit-est -i uniq.F.fasta -o xxx -c $CDHIT -T 0 -M 0 -g 1 -d 100 &>cdhit.log
	mawk '{if ($1 ~ /Cl/) clus = clus + 1; else  print $3 "\t" clus}' xxx.clstr | sed 's/[>dDocent_Contig_,...]//g' | sort -g -k1 > sort.contig.cluster.ids
	paste sort.contig.cluster.ids totaluniqseq > contig.cluster.totaluniqseq
	sort -k2,2 -g contig.cluster.totaluniqseq | sed -e 's/NNNNNNNNNN/	/g' > rcluster
	#CD-hit output is converted to rainbow format
	rainbow div -i rcluster -o rbdiv.out -f 0.5 -K 10
	if [ "$ATYPE" == "PE" ]; then
		rainbow merge -o rbasm.out -a -i rbdiv.out -r 2 -N10000 -R10000 -l 20 -f 0.75
	else
		rainbow merge -o rbasm.out -a -i rbdiv.out
	fi
	#This AWK code replaces rainbow's contig selection perl script
	cat rbasm.out <(echo "E") |sed 's/[0-9]*:[0-9]*://g' | mawk ' {
		if (NR == 1) e=$2;
		else if ($1 ~/E/ && lenp > len1) {c=c+1; print ">dDocent_Contig_" e "\n" seq2 "NNNNNNNNNN" seq1; seq1=0; seq2=0;lenp=0;e=$2;fclus=0;len1=0;freqp=0;lenf=0}
		else if ($1 ~/E/ && lenp <= len1) {c=c+1; print ">dDocent_Contig_" e "\n" seq1; seq1=0; seq2=0;lenp=0;e=$2;fclus=0;len1=0;freqp=0;lenf=0}
		else if ($1 ~/C/) clus=$2;
		else if ($1 ~/L/) len=$2;
		else if ($1 ~/S/) seq=$2;
		else if ($1 ~/N/) freq=$2;
		else if ($1 ~/R/ && $0 ~/0/ && $0 !~/1/ && len > lenf) {seq1 = seq; fclus=clus;lenf=len}
		else if ($1 ~/R/ && $0 ~/0/ && $0 ~/1/) {seq1 = seq; fclus=clus; len1=len}
		else if ($1 ~/R/ && $0 ~!/0/ && freq > freqp && len >= lenp || $1 ~/R/ && $0 ~!/0/ && freq == freqp && len > lenp) {seq2 = seq; lenp = len; freqp=freq}
        }' > rainbow.fasta

	seqtk seq -r rainbow.fasta > rainbow.RC.fasta
	mv rainbow.RC.fasta rainbow.fasta

	#The rainbow assembly is checked for overlap between newly assembled Forward and Reverse reads using the software PEAR
	sed -e 's/NNNNNNNNNN/	/g' rainbow.fasta | cut -f1 | gawk 'BEGIN {RS = ">" ; FS = "\n"} NR > 1 {print "@"$1"\n"$2"\n+""\n"gensub(/./, "I", "g", $2)}' > ref.F.fq
	sed -e 's/NNNNNNNNNN/	/g' rainbow.fasta | cut -f2 | gawk 'BEGIN {RS = ">" ; FS = "\n"} NR > 1 {print "@"$1"\n"$2"\n+""\n"gensub(/./, "I", "g", $2)}' > ref.R.fq

	seqtk seq -r ref.R.fq > ref.RC.fq
	mv ref.RC.fq ref.R.fq
	LENGTH=$(mawk '!/>/' rainbow.fasta | mawk '(NR==1||length<shortest){shortest=length} END {print shortest}')
	LENGTH=$(( $LENGTH * 5 / 4))
	

	pearRM -f ref.F.fq -r ref.R.fq -o overlap -p 0.001 -j $NUMProc -n $LENGTH

	rm ref.F.fq ref.R.fq

	mawk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' overlap.assembled.fastq > overlap.fasta
	mawk '/>/' overlap.fasta > overlap.loci.names
	mawk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' overlap.unassembled.forward.fastq > other.F
	mawk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' overlap.unassembled.reverse.fastq > other.R
	paste other.F other.R | mawk '{if ($1 ~ />/) print $1; else print $0}' | sed 's/	/NNNNNNNNNN/g' > other.FR

	cat other.FR overlap.fasta > totalover.fasta

	rm *.F *.R
fi
if [[ "$ATYPE" != "PE" && "$ATYPE" != "RPE" ]]; then
	cp uniq.fasta totalover.fasta
fi
cd-hit-est -i totalover.fasta -o reference.fasta.original -M 0 -T 0 -c $simC

sed -e 's/^C/NC/g' -e 's/^A/NA/g' -e 's/^G/NG/g' -e 's/^T/NT/g' -e 's/T$/TN/g' -e 's/A$/AN/g' -e 's/C$/CN/g' -e 's/G$/GN/g' reference.fasta.original > reference.fasta


samtools faidx reference.fasta
bwa index reference.fasta

}

##Create alignment intervals
##This takes advantage of the fact that RAD loci are very discrete.  Instead of calculating intervals for every BAM file,
##this function merges all BAM files together.  This overall BAM file 
##is used to create a single list of intervals, saving a large amount of computational time.

CreateIntervals()
{
samtools merge -@$NUMProc -b bamlist.list -f cat-RRG.bam &>/dev/null
samtools index cat-RRG.bam 
wait
bamToBed -i cat-RRG.bam | bedtools merge -i - > mapped.bed
}

#This checks that dDocent has detected the proper number of individuals and exits if incorrect
GetInfo(){
echo "$NumInd individuals are detected. Is this correct? Enter yes or no and press [ENTER]"

read Indcorrect

if [ "$Indcorrect" == "no" ]; then
        echo "Please double check that all fastq files are named Ind01.F.fq.gz and Ind01.R.fq.gz"
        exit 1
elif [ "$Indcorrect" == "yes" ]; then
            echo "Proceeding with $NumInd individuals"
else
        echo "Incorrect Input"
        exit 1
fi

#Tries to get number of processors, if not asks user
NUMProc=( `grep -c ^processor /proc/cpuinfo 2> /dev/null` ) 
NUMProc=$(($NUMProc + 0)) 

echo "dDocent detects $NUMProc processors available on this system."
echo "Please enter the maximum number of processors to use for this analysis."
        read NUMProc
        
if [ $NUMProc -lt 1 ]; then
        echo "Incorrect. Please enter the number of processing cores on this computer"
        read NUMProc
fi                
if [ $NUMProc -lt 1 ]; then
        echo "Incorrect input, exiting"
        exit 1
fi

#Tries to get maximum system memory, if not asks user
MAXMemory=$(($(grep -Po '(?<=^MemTotal:)\s*[0-9]+' /proc/meminfo | tr -d " ") / 1048576))G

echo "dDocent detects $MAXMemory maximum memory available on this system."
echo "Please enter the maximum memory to use for this analysis. The size can be postfixed with 
K, M, G, T, P, k, m, g, t, or p which would multiply the size with 1024, 1048576, 1073741824, 
1099511627776, 1125899906842624, 1000, 1000000, 1000000000, 1000000000000, or 1000000000000000 respectively."
echo "For example, to limit dDocent to ten gigabytes, enter 10G or 10g"
        read MAXMemory

while [[ -z $MAXMemory ]];
	do
	echo "Incorrect input"
	echo -e "Please enter the maximum memory to use for this analysis. The size can be postfixed with K, M, G, T, P, k, m, g, t, or p which would multiply the size with 1024, 1048576, 1073741824, 1099511627776, 1125899906842624, 1000, 1000000, 1000000000, 1000000000000, or 1000000000000000 respectively."
	echo -e "This option does not work with all distributions of Linux.  If runs are hanging at variant calling, enter 0"
	echo -e "Then press [ENTER]"
	read MAXMemory
	done

#Asks if user wants to trim reads.  This allows this part of the pipeline to be skipped during subsequent analyses
echo -e "\nDo you want to quality trim your reads?" 
echo "Type yes or no and press [ENTER]?"

read TRIM

#Asks if user wants to perform an assembly.  This allows this part of the pipeline to be skipped during subsequent analyses

echo -e "\nDo you want to perform an assembly?"
echo "Type yes or no and press [ENTER]."

read ASSEMBLY

if [ "$ASSEMBLY" == "no" ]; then
        echo -e "\nReference contigs need to be in a file named reference.fasta\n"
        sleep 1
else
	echo -e "What type of assembly would you like to perform?  Enter SE for single end, PE for paired-end, RPE for paired-end sequencing for RAD protocols with random shearing, or OL for paired-end sequencing that has substantial overlap."
	echo -e "Then press [ENTER]"
	read ATYPE

	while [[ $ATYPE != "SE" && $ATYPE != "PE" && $ATYPE != "OL" && $ATYPE != "RPE" ]];
	do
	echo "Incorrect input"
	echo -e "What type of assembly would you like to perform?  Enter SE for single end, PE for paired-end, RPE for paired-end sequencing for RAD protocols with random shearing, or OL for paired-end sequencing that has substantial overlap."
	echo -e "Then press [ENTER]"
	read ATYPE
	done
fi
#If performing de novo assembly, asks if the user wants to enter a different -c value
if [ "$ASSEMBLY" == "yes" ]; then
    echo "Reads will be assembled with Rainbow"
    echo "CD-HIT will cluster reference sequences by similarity. The -c parameter (% similarity to cluster) may need to be changed for your taxa."
    echo "Would you like to enter a new c parameter now? Type yes or no and press [ENTER]"
    read optC
    if [ "$optC" == "no" ]; then
            echo "Proceeding with default 0.9 value."
            simC=0.9
        elif [ "$optC" == "yes" ]; then
            echo "Please enter new value for c. Enter in decimal form (For 90%, enter 0.9)"
            read newC
            simC=$newC
        else
            echo "Incorrect input. Proceeding with the default value."
            simC=0.9
        fi
fi

#Asks if user wants to map reads and change default mapping variables for BWA
echo "Do you want to map reads?  Type yes or no and press [ENTER]"
read MAP
if [ "$MAP" == "no" ]; then
        echo "Mapping will not be performed"
        optA=1
    	optB=4
    	optO=6
        else
                echo "BWA will be used to map reads.  You may need to adjust -A -B and -O parameters for your taxa."
                echo "Would you like to enter a new parameters now? Type yes or no and press [ENTER]"
                read optq

        if [ "$optq" == "yes" ]; then
        echo "Please enter new value for A (match score).  It should be an integer.  Default is 1."
        read newA
        optA=$newA
                echo "Please enter new value for B (mismatch score).  It should be an integer.  Default is 4."
        read newB
        optB=$newB
                echo "Please enter new value for O (gap penalty).  It should be an integer.  Default is 6."
        read newO
        optO=$newO
        else
                echo "Proceeding with default values for BWA read mapping."
                optA=1
                optB=4
                optO=6
        fi
fi

#Does user wish to call SNPs?
echo "Do you want to use FreeBayes to call SNPs?  Please type yes or no and press [ENTER]"
read SNP

while [[ $SNP != "yes" && $SNP != "no" ]];
	do
	echo "Incorrect input"
	echo -e "Do you want to use FreeBayes to call SNPs?  Please type yes or no and press [ENTER]"
	read SNP
	done

#Asks user for email address to notify when analysis is complete
echo ""
echo "Please enter your email address.  dDocent will email you when it is finished running."
echo "Don't worry; dDocent has no financial need to sell your email address to spammers."
read MAIL
echo ""
echo ""

if [ "$ASSEMBLY" == "no" ]; then
#Prints instructions on how to move analysis to background and disown process
echo "At this point, all configuration information has been entered and dDocent may take several hours to run." 
echo "It is recommended that you move this script to a background operation and disable terminal input and output."
echo "All data and logfiles will still be recorded."
echo "To do this:"
echo "Press control and Z simultaneously"
echo "Type 'bg' without the quotes and press enter"
echo "Type 'disown -h' again without the quotes and press enter"
echo ""
echo "Now sit back, relax, and wait for your analysis to finish."
fi

if [ "$ASSEMBLY" == "yes" ]; then
echo "dDocent will require input during the assembly stage.  Please wait until prompt says it is safe to move program to the background."
fi
}

#Actually starts program
if [ -n "$1" ]; then
	main $1 2>&1 | tee -a dDocent_main.LOG #Log all output
else
	main 2>&1 | tee -a dDocent_main.LOG  #Log all output
fi

#Compress Large Leftover files
gzip -f concat.fasta concat.seq rcluster rbdiv.out rbasm.out rainbow.fasta reference.fasta.original uniq.seqs uniq.fasta totaluniqseq uniq.F.fasta uniq.RC.fasta 2> /dev/null &

