#!/bin/bash

#########dDocent 1.0##################################################################################


#This script serves as a bash wrapper to QC, assemble, map, and call SNPs from double digest RAD data.
#It requires that your raw data are split up by tagged individual and follow the naming convenction of:

#Sample1.F.fq and Sample1.R.fq

#This is the non-interactive version, you must enter 5 variables with the shell script.

#Correct usage is: sh NI.dDocent.sh [c] [A] [B] [O] [email address]

#All five parameters must be specificed

#The first parameter is the c parameter:

#This paramter changes the similarity threshhold to cluster sequences.  
#This is sensitive to the amount of natural variation in your population.
#Set it too high and you will have homologous loci that will stay seperate in your analysis and most likely be lost.
#Set it too low and you will possibly cluster paralogs 

#The other parameters control mapping

#-A controls the match score, -B controls the mismatch score, and -O controls the Gap Penalty
#Optimal mapping will likely involved tweaking these parameters, especially -0

#I have commonly used -A 2 -B 4 -O 2, -A 4 -B 7 -O 4, -A 3, -B 5, -O 4 for examples

###################################################################################################

if [ ! $# == 5 ]; then
	echo "You must enter all 5 parameters."
	echo -e "Correct usage is: sh NI.dDocent.sh [c] [A] [B] [O] [email address] \nAll five parameters must be specificed\n\nThe first parameter is the c parameter: \nThis paramter changes the similarity threshhold to cluster sequences. \nThis is sensitive to the amount of natural variation in your population. \nSet it too high and you will have homologous loci that will stay seperate in your analysis and most likely be lost.\nSet it too low and you will possibly cluster paralogs. \n" 
	echo -e "The next 3  parameters control mapping: \n-A controls the match score, -B controls the mismatch score, and -O controls the Gap Penalty\nOptimal mapping will likely involved tweaking these parameters, especially -O \n"
	echo -e "I have commonly used -A 2 -B 4 -O 2, -A 4 -B 7 -O 4, -A 3, -B 5, -O 4 \n"
	echo -e "The last parameter is your email address.  dDocent will email you when your analysis is complete. \n"
	echo -e "Default usage would be sh NI.dDocent.sh 0.9 1 4 6 jpuritz@gmail.com"
	exit 1
fi

NumInd=$(ls *.F.fq | wc -l)
NumInd=$(($NumInd - 0))

echo -e "dDocent 1.0 (Non-interactive Version) by J. Puritz for Gold lab \n"
echo -e "Contact jpuritz@gmail.com with any problems \n\n "
if [ $NumInd -gt 9 ]
        then
	MinAll=0.05
	MaxSize=9
        else
	MinAll=$(echo "scale=2; 1 / (2 * $NumInd) " | bc)
	MaxSize=$(( $NumInd - 1 ))
fi


ls *.F.fq > namelist
sed -i 's/.F.fq//g' namelist
NAMES=( `cat "namelist" `)

ls -S *.F.fq > sizelist
sed -i 's/.F.fq//g' sizelist
SIZE=( `cat "sizelist" `)

	echo "Removing the _1 character and replacing with /1 in the name of every sequence"
	for i in "${NAMES[@]}"
	do	
	sed -e 's:_2$:/2:g' $i.R.fq > $i.Ra.fq
	sed -e 's:_1$:/1:g' $i.F.fq > $i.Fa.fq
	mv $i.Ra.fq $i.R.fq
	mv $i.Fa.fq $i.F.fq
	done

TrimReads () 
{ for i in "${NAMES[@]}"
do
echo "Trimming $i"
trim_galore --paired -q 10 --length 20 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATATCGTATGCCGTCTTCTGCTTG -a2 GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCG --stringency 10 $i.F.fq $i.R.fq 2> $i.trim.log
mv $i.F_val_1.fq $i.R1.fq
mv $i.R_val_2.fq $i.R2.fq
done
}

#Use Rainbow to cluster and assemble reads

setupRainbow ()
{ echo "Concatenating F and R reads of up to 10 individuals for assembly"
cat ${SIZE[0]}.F.fq > forward
cat ${SIZE[0]}.R.fq > reverse

for ((i = 1; i <= $MaxSize; i++));
do
cat ${SIZE[$i]}.F.fq >> forward
cat ${SIZE[$i]}.R.fq >> reverse
done

seqtk seq -r forward > forwardRC
mergefq.pl reverse forwardRC concat.fq

#Use Rainbow to cluster and assemble reads
echo "Using rainbow to cluster and assemble reads"
rainbow cluster -1 concat.fq -m 6 > cat.rbcluster.out 2> log
rainbow div -i cat.rbcluster.out -o cat.rbdiv.out -f $MinAll
rainbow merge -a -i cat.rbdiv.out -o cat.rbasm.out
select_best_rbcontig.pl cat.rbasm.out > rainbow
cat rainbow | sed s/./N/96 | sed s/./N/97 | sed s/./N/98 | sed s/./N/99 | sed s/./N/100 | sed s/./N/101 | sed s/./N/102 | sed s/./N/103 | sed s/./N/104 | sed s/./N/105 > rainbowN
echo "Now using cd-hit to cluster reference sequences by similarity. The -c parameter (% similarity to cluster) may need to be changed for your taxa"
cd-hit-est -i rainbowN -o referencegenome -T 0 -c $1 -M 0 -l 30 &>cdhit.log
}

	echo "Trimming reads and simultaneously assemblying reference sequences"	
	TrimReads & 2> trim.log; setupRainbow $1 2> rainbow.log
	wait


##Use BWA to map reads to assembly
echo "Using BWA to map reads.  You may need to adjust -A -B and -O parameters for your taxa."
bwa0.7 index -a bwtsw referencegenome &> index.log

for i in "${NAMES[@]}"
do
bwa0.7 mem referencegenome $i.R1.fq $i.R2.fq -t 32 -a -T 10 -A $2 -B $3 -O $4 > $i.sam 2> bwa.$i.log
done

##Convert Sam to Bam and remove low quality, ambiguous mapping
for i in "${NAMES[@]}"
do
samtools view -bT referencegenome -q1 $i.sam > $i.bam 2>$i.bam.log
samtools sort $i.bam $i
done

samtools faidx referencegenome

#Calling of SNPs from two samples must have a minimum read of depth of 10 and below 200 with a minimum quality score of 20
echo "Using samtools to pileup reads"
samtools mpileup -D -f referencegenome *.bam >mpileup 2> mpileup.log
echo "Using VarScan2 to call SNPs with at least 5 reads (within 1 individual), 95% probability, and at least 2 reads for the minor allele"
java -jar /usr/local/bin/VarScan.v2.3.5.jar mpileup2snp mpileup --output-vcf --min-coverage 5 --strand-filter 0 --min-var-freq 0.1 --p-value 0.05 >SNPS.vcf 2>varscan.log
echo "Using VCFtools to parse SNPS.vcf for SNPS that are not indels and are called in at least 90% of individuals"
vcftools --vcf SNPS.vcf --geno 0.9 --out Final --counts --recode --non-ref-af 0.001 --remove-indels &>VCFtools.log

tail Final.log	

if [ ! -d "logfiles" ];	then
mkdir logfiles
fi

mv *.log *.txt log ./logfiles
echo `pwd` `date` | mailx -s "It's Finished" $5
