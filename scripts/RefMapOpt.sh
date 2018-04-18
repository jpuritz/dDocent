#!/usr/bin/env bash
export LC_ALL=en_US.UTF-8

if [[ -z "$7" ]]; then
echo "Usage is RefMapOpt minK1 maxK1 minK2 maxK2 cluster_similarity Num_of_Processors Assembly_Type optional_list_of_individuals"
exit 1
fi

if ! awk --version | fgrep -v GNU &>/dev/null; then
         awk=gawk
    else
         awk=awk
fi

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

ATYPE=$7
if [[ $ATYPE != "SE" && $ATYPE != "PE" && $ATYPE != "OL" && $ATYPE != "HYB" && $ATYPE != "ROL" ]]; then
echo "Usage is RefMapTest minK1 maxK1 minK2 maxK2 cluster_similarity Num_of_Processors Assembly_Type optional_list_of_individuals"
echo "Please make sure to choose assembly type."
exit 1
fi

NUMProc=$6
ls *.F.fq.gz > namelist
sed -i'' -e 's/.F.fq.gz//g' namelist
NAMES=( `cat "namelist" `)

#This code checks for trimmed sequence files
TEST=$(ls *.R1.fq.gz 2> /dev/null | wc -l )
if [ "$TEST" -gt 0 ]; then
	echo "Trimmed sequences found, proceeding with optimization."
else
	echo -e "\nRefMapOpt.sh requires that you have trimmed sequence files.\nPlease include trimmed sequence files with the .R1.fq.gz and .R2.fq.gz naming convention."
	echo "dDocent will create these for you"
	echo "Please rerun RefMapOpt.sh after trimming sequence files"
exit 1
fi


Reference(){

AWK1='BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}'
AWK2='!/>/'
AWK3='!/NNN/'
PERLT='while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}'
SED1='s/^[ \t]*//'
SED2='s/\s/\t/g'
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
	if [[ "$ATYPE" == "OL" || "$ATYPE" == "ROL" ]]; then
	#If OL assembly, dDocent assumes that the marjority of PE reads will overlap, so the software PEAR is used to merge paired reads into single reads
		for i in "${NAMES[@]}";
        		do
        		zcat $i.R.fq.gz | head -2 | tail -1 >> lengths.txt
        		done	
        MaxLen=$(mawk '{ print length() | "sort -rn" }' lengths.txt| head -1)
		LENGTH=$(( $MaxLen / 3))
		for i in "${NAMES[@]}"
			do
			echo "Running PEAR on sample" $i
			pearRM -f $i.F.fq.gz -r $i.R.fq.gz -o $i -j $NUMProc -n $LENGTH &> $i.kopt.log
			done
		cat namelist | parallel --no-notice -j $NUMProc "mawk '$AWK1' {}.assembled.fastq | mawk '$AWK2' | perl -e '$PERLT' > {}.uniq.seqs"
	fi
	if [ "$ATYPE" == "HYB" ]; then
	#If HYB assembly, dDocent assumes some PE reads will overlap but that some will not, so the OL method performed and remaining reads are then put through PE method
		for i in "${NAMES[@]}";
        		do
        		zcat $i.R.fq.gz | head -2 | tail -1 >> lengths.txt
        		done	
        	MaxLen=$(mawk '{ print length() | "sort -rn" }' lengths.txt| head -1)
		LENGTH=$(( $MaxLen / 3))
		for i in "${NAMES[@]}"
			do
			pearRM -f $i.F.fq.gz -r $i.R.fq.gz -o $i -j $NUMProc -n $LENGTH &>kopt.log
			done
		cat namelist | parallel --no-notice -j $NUMProc "mawk '$AWK1' {}.assembled.fastq | mawk '$AWK2' | perl -e '$PERLT' > {}.uniq.seqs"
		
		cat namelist | parallel --no-notice -j $NUMProc "cat {}.unassembled.forward.fastq | mawk '$AWK1' | mawk '$AWK2' > {}.forward"
		cat namelist | parallel --no-notice -j $NUMProc "cat {}.unassembled.reverse.fastq | mawk '$AWK1' | mawk '$AWK2' > {}.reverse"
		cat namelist | parallel --no-notice -j $NUMProc "paste -d '-' {}.forward {}.reverse | mawk '$AWK3'| sed 's/-/NNNNNNNNNN/' | perl -e '$PERLT' > {}.uniq.ua.seqs"
		rm *.forward
		rm *.reverse
	fi		
	
fi

#Create a data file with the number of unique sequences and the number of occurrences
if [[ "$ATYPE" == "RPE" || "$ATYPE" == "ROL" ]]; then
	parallel --no-notice mawk -v x=$1 \''$1 >= x'\' ::: *.uniq.seqs | cut -f2 |  sort | uniq -c -w $FRL | sed -e 's/^[ \t]*//' | sed -e 's/\s/\t/g' | mawk -v x=$2 '$1 >= x' > uniq.k.$1.c.$2.seqs
else
	parallel --no-notice mawk -v x=$1 \''$1 >= x'\' ::: *.uniq.seqs | cut -f2 | perl -e 'while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}' | mawk -v x=$2 '$1 >= x' > uniq.k.$1.c.$2.seqs
fi
cut -f2 uniq.k.$1.c.$2.seqs > totaluniqseq
mawk '{c= c + 1; print ">dDocent_Contig_" c "\n" $1}' totaluniqseq > uniq.full.fasta
LENGTH=$(mawk '!/>/' uniq.full.fasta  | mawk '(NR==1||length<shortest){shortest=length} END {print shortest}')
LENGTH=$(($LENGTH * 3 / 4))
$awk 'BEGIN {RS = ">" ; FS = "\n"} NR > 1 {print "@"$1"\n"$2"\n+""\n"gensub(/./, "I", "g", $2)}' uniq.full.fasta > uniq.fq
java -jar $TRIMMOMATIC SE -threads $NUMProc -phred33 uniq.fq uniq.fq1 ILLUMINACLIP:$ADAPTERS:2:30:10 MINLEN:$LENGTH &>/dev/null
mawk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' uniq.fq1 > uniq.fasta
mawk '!/>/' uniq.fasta > totaluniqseq
rm uniq.fq*

#If this is a PE assebmle
if [[ "$ATYPE" == "PE" || "$ATYPE" == "RPE" ]]; then
	#Reads are first clustered using only the Forward reads using CD-hit instead of rainbow
	sed -e 's/NNNNNNNNNN/	/g' uniq.fasta | cut -f1 > uniq.F.fasta
	CDHIT=$(python -c "print (max("$3" - 0.1,0.8))")
	cd-hit-est -i uniq.F.fasta -o xxx -c $CDHIT -T 0 -M 0 -g 1 -d 100 &>cdhit.log
	mawk '{if ($1 ~ /Cl/) clus = clus + 1; else  print $3 "\t" clus}' xxx.clstr | sed 's/[>dDocent_Contig_,...]//g' | sort -g -k1 > sort.contig.cluster.ids
	paste sort.contig.cluster.ids totaluniqseq > contig.cluster.totaluniqseq
	sort -k2,2 -g contig.cluster.totaluniqseq | sed -e 's/NNNNNNNNNN/	/g' > rcluster
	#CD-hit output is converted to rainbow format
	rainbow div -i rcluster -o rbdiv.out -f 0.5 -K 10
	rainbow merge -o rbasm.out -a -i rbdiv.out -r 2 -N10000 -R10000 -l 20 -f 0.75
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
	sed -e 's/NNNNNNNNNN/	/g' rainbow.fasta | cut -f1 | $awk 'BEGIN {RS = ">" ; FS = "\n"} NR > 1 {print "@"$1"\n"$2"\n+""\n"gensub(/./, "I", "g", $2)}' > ref.F.fq
	sed -e 's/NNNNNNNNNN/	/g' rainbow.fasta | cut -f2 | $awk 'BEGIN {RS = ">" ; FS = "\n"} NR > 1 {print "@"$1"\n"$2"\n+""\n"gensub(/./, "I", "g", $2)}' > ref.R.fq

	seqtk seq -r ref.R.fq > ref.RC.fq
	mv ref.RC.fq ref.R.fq
	LENGTH=$(mawk '!/>/' rainbow.fasta | mawk '(NR==1||length<shortest){shortest=length} END {print shortest}')
	LENGTH=$(( $LENGTH * 5 / 4))
	
	pearRM -f ref.F.fq -r ref.R.fq -o overlap -p 0.001 -j $NUMProc -n $LENGTH &>kopt.log

	rm ref.F.fq ref.R.fq

	mawk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' overlap.assembled.fastq > overlap.fasta
	mawk '/>/' overlap.fasta > overlap.loci.names
	mawk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' overlap.unassembled.forward.fastq > other.F
	mawk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' overlap.unassembled.reverse.fastq > other.R
	paste other.F other.R | mawk '{if ($1 ~ />/) print $1; else print $0}' | sed 's/	/NNNNNNNNNN/g' > other.FR

	cat other.FR overlap.fasta > totalover.fasta

	rm *.F *.R
fi

if [[ "$ATYPE" == "HYB" ]];then
	parallel --no-notice mawk -v x=$1 \''$1 >= x'\' ::: *.uniq.ua.seqs | cut -f2 | perl -e 'while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}' | mawk -v x=$2 '$1 >= x' > uniq.k.$1.c.$2.ua.seqs
	AS=$(cat uniq.k.$1.c.$2.ua.seqs | wc -l)
	if [ "$AS" -gt 1 ]; then
		cut -f2 uniq.k.$1.c.$2.ua.seqs > totaluniqseq.ua
		mawk '{c= c + 1; print ">dDocent_Contig_" c "\n" $1}' totaluniqseq.ua > uniq.full.ua.fasta
		LENGTH=$(mawk '!/>/' uniq.full.ua.fasta  | mawk '(NR==1||length<shortest){shortest=length} END {print shortest}')
		LENGTH=$(($LENGTH * 3 / 4))
		$awk 'BEGIN {RS = ">" ; FS = "\n"} NR > 1 {print "@"$1"\n"$2"\n+""\n"gensub(/./, "I", "g", $2)}' uniq.full.ua.fasta > uniq.ua.fq
		java -jar $TRIMMOMATIC SE -threads $NUMProc -phred33 uniq.ua.fq uniq.ua.fq1 ILLUMINACLIP:$ADAPTERS:2:30:10 MINLEN:$LENGTH &>/dev/null
		mawk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' uniq.ua.fq1 > uniq.ua.fasta
		mawk '!/>/' uniq.ua.fasta > totaluniqseq.ua
		rm uniq.ua.fq*
		#Reads are first clustered using only the Forward reads using CD-hit instead of rainbow
		sed -e 's/NNNNNNNNNN/	/g' uniq.ua.fasta | cut -f1 > uniq.F.ua.fasta
		CDHIT=$(python -c "print(max("$3" - 0.1,0.8))")
		cd-hit-est -i uniq.F.ua.fasta -o xxx -c $CDHIT -T 0 -M 0 -g 1 -d 100 &>cdhit.log
		mawk '{if ($1 ~ /Cl/) clus = clus + 1; else  print $3 "\t" clus}' xxx.clstr | sed 's/[>dDocent_Contig_,...]//g' | sort -g -k1 > sort.contig.cluster.ids.ua
		paste sort.contig.cluster.ids.ua totaluniqseq.ua > contig.cluster.totaluniqseq.ua
		sort -k2,2 -g contig.cluster.totaluniqseq.ua | sed -e 's/NNNNNNNNNN/	/g' > rcluster.ua
		#CD-hit output is converted to rainbow format
		rainbow div -i rcluster.ua -o rbdiv.ua.out -f 0.5 -K 10
		rainbow merge -o rbasm.ua.out -a -i rbdiv.ua.out -r 2 -N10000 -R10000 -l 20 -f 0.75
		#This AWK code replaces rainbow's contig selection perl script
		cat rbasm.ua.out <(echo "E") |sed 's/[0-9]*:[0-9]*://g' | mawk ' {
			if (NR == 1) e=$2;
			else if ($1 ~/E/ && lenp > len1) {c=c+1; print ">dDocent_Contig_UA_" e "\n" seq2 "NNNNNNNNNN" seq1; seq1=0; seq2=0;lenp=0;e=$2;fclus=0;len1=0;freqp=0;lenf=0}
			else if ($1 ~/E/ && lenp <= len1) {c=c+1; print ">dDocent_Contig_UA_" e "\n" seq1; seq1=0; seq2=0;lenp=0;e=$2;fclus=0;len1=0;freqp=0;lenf=0}
			else if ($1 ~/C/) clus=$2;
			else if ($1 ~/L/) len=$2;
			else if ($1 ~/S/) seq=$2;
			else if ($1 ~/N/) freq=$2;
			else if ($1 ~/R/ && $0 ~/0/ && $0 !~/1/ && len > lenf) {seq1 = seq; fclus=clus;lenf=len}
			else if ($1 ~/R/ && $0 ~/0/ && $0 ~/1/) {seq1 = seq; fclus=clus; len1=len}
			else if ($1 ~/R/ && $0 ~!/0/ && freq > freqp && len >= lenp || $1 ~/R/ && $0 ~!/0/ && freq == freqp && len > lenp) {seq2 = seq; lenp = len; freqp=freq}
			}' > rainbow.ua.fasta
	
		seqtk seq -r rainbow.ua.fasta > rainbow.RC.fasta
		mv rainbow.RC.fasta rainbow.ua.fasta
	
		cat rainbow.ua.fasta uniq.fasta > totalover.fasta

	fi
fi

if [[ "$ATYPE" != "PE" && "$ATYPE" != "RPE" && "$ATYPE" != "HYB" ]]; then
	cp uniq.fasta totalover.fasta
fi

cd-hit-est -i totalover.fasta -o reference.fasta.original -M 0 -T 0 -c $3 &>cdhit.log

sed -e 's/^C/NC/g' -e 's/^A/NA/g' -e 's/^G/NG/g' -e 's/^T/NT/g' -e 's/T$/TN/g' -e 's/A$/AN/g' -e 's/C$/CN/g' -e 's/G$/GN/g' reference.fasta.original > reference.fasta


samtools faidx reference.fasta
bwa index reference.fasta
}


ls -S *.F.fq.gz > namelist
sed -i'' -e 's/.F.fq.gz//g' namelist

if [[ -z "$8" ]]; then
NUMIND=$(cat namelist | wc -l)
NUM1=$(($NUMIND / 10))
NUM2=$(($NUMIND - $NUM1))
NUM3=$(($NUM2 - $NUM1))
cat namelist | head -$NUM2 | tail -$NUM3 > newlist
NAMESR=( `cat "newlist" `)
LEN=${#NAMESR[@]}
LEN=$(($LEN - 1))
echo "5" >randlist
	for ((rr = 1; rr<=50; rr++));
	do
	INDEX=$[ 1 + $[ RANDOM % $LEN ]]
		if grep -q ${NAMESR[$INDEX]} randlist;
		then x=x
		else
	echo ${NAMESR[$INDEX]} >> randlist
		fi
	done

RANDNAMES=( `mawk '!/^5$/' "randlist" | head -20 `)
else
RANDNAMES=( `cat "$8" `)
fi

echo -e "Cov\tNon0Cov\tContigs\tMeanContigsMapped\tK1\tK2\tSUM Mapped\tSUM Properly\tMean Mapped\tMean Properly\tMisMatched" > mapping.results

for ((r = $1; r <= $2; r++));
do
	for ((j = $3; j <= $4; j++));
	do
	PP=$(($r + $j))
	if [ "$PP" != "2" ]; then
	Reference $r $j $5
    	rm lengths.txt &> /dev/null
    		for k in "${RANDNAMES[@]}";
    		do
    		if [ -f "$k.R.fq.gz" ]; then
    		zcat $k.R.fq.gz | head -2 | tail -1 >> lengths.txt
    		fi
    		done	
    	if [ -f "lengths.txt" ]; then
    	MaxLen=$(mawk '{ print length() | "sort -rn" }' lengths.txt| head -1)
    	INSERT=$(($MaxLen * 2 ))
    	INSERTH=$(($INSERT + 100 ))
    	INSERTL=$(($INSERT - 100 ))
    	SD=$(($INSERT / 10))
		fi
    	#BWA for mapping for all samples
    	rm $r.$j.results 2>/dev/null
    		for k in "${RANDNAMES[@]}"
    		do
		if [[ "$ATYPE" == "OL" || "$ATYPE" == "HYB"  || "$ATYPE" == "ROL" ]]; then
			if [ -f "$k.R2.fq.gz" ]; then
				bwa mem reference.fasta $k.R1.fq.gz $k.R2.fq.gz -L 20,5 -t 32 -a -M -T 10 -A1 -B 3 -O 5 -R "@RG\tID:$k\tSM:$k\tPL:Illumina" 2> bwa.$k.log | mawk '!/\t[2-9].[SH].*/' | mawk '!/[2-9].[SH]\t/' | samtools view -@32 -q 1 -SbT reference.fasta - > $k.bam
			else
				bwa mem reference.fasta $k.R1.fq.gz -L 20,5 -t 32 -a -M -T 10 -A1 -B 3 -O 5 -R "@RG\tID:$k\tSM:$k\tPL:Illumina" 2> bwa.$k.log | mawk '!/\t[2-9].[SH].*/' | mawk '!/[2-9].[SH]\t/' | samtools view -@32 -q 1 -SbT reference.fasta - > $k.bam
			fi
		else
			if [ -f "$k.R2.fq.gz" ]; then
    			bwa mem reference.fasta $k.R1.fq.gz $k.R2.fq.gz -L 20,5 -I $INSERT,$SD,$INSERTH,$INSERTL -t 32 -a -M -T 10 -A 1 -B 3 -O 5 -R "@RG\tID:$k\tSM:$k\tPL:Illumina" 2> bwa.$k.log | mawk '!/\t[2-9].[SH].*/' | mawk '!/[2-9].[SH]\t/' | samtools view -@32 -q 1 -SbT reference.fasta - > $k.bam
    		else
    			bwa mem reference.fasta $k.R1.fq.gz -L 20,5 -t 32 -a -M -T 10 -A 1 -B 3 -O 5 -R "@RG\tID:$k\tSM:$k\tPL:Illumina" 2> bwa.$k.log | mawk '!/\t[2-9].[SH].*/' | mawk '!/[2-9].[SH]\t/' | samtools view -@32 -q 1 -SbT reference.fasta - > $k.bam
    		fi
    	fi
		samtools sort -@24 $k.bam -o $k.bam 
		samtools index $k.bam
    		MM=$(samtools flagstat $k.bam | grep -E 'mapped \(|properly' | cut -f1 -d '+' | tr -d '\n')
    		CM=$(samtools idxstats $k.bam | mawk '$3 > 0' | wc -l)
    		CC=$(samtools idxstats $k.bam | mawk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }')
		DD=$(samtools idxstats $k.bam | mawk '$3 >0' | mawk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }')
		BM=$(samtools flagstat $k.bam | grep mapQ | cut -f1 -d ' ')
		echo -e "$MM\t$CM\t$CC\t$DD\t$BM" >> $r.$j.results
    		done
    	SUMM=$(mawk '{ sum+=$1} END {print sum}' $r.$j.results)
    	SUMPM=$(mawk '{ sum+=$2} END {print sum}' $r.$j.results)
    	AVEM=$(mawk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' $r.$j.results)
    	AVEP=$(mawk '{ sum += $2; n++ } END { if (n > 0) print sum / n; }' $r.$j.results)
    	AC=$(mawk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' $r.$j.results)
    	CON=$(mawk '/>/' reference.fasta | wc -l)
	CCC=$(mawk '{ sum += $4; n++ } END { if (n > 0) print sum / n; }' $r.$j.results)
	DDD=$(mawk '{ sum += $5; n++ } END { if (n > 0) print sum / n; }' $r.$j.results)
	BBM=$(mawk '{ sum += $6; n++ } END { if (n > 0) print sum / n; }' $r.$j.results)
    	echo -e "$CCC\t$DDD\t$CON\t$AC\t$r\t$j\t$SUMM\t$SUMPM\t$AVEM\t$AVEP\t$BBM" >> mapping.results
	fi
	done
    		
done
