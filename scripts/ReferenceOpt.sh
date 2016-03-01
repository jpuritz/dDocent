#!/usr/bin/env bash

if [[ -z "$6" ]]; then
echo "Usage is sh ReferenceOpt.sh minK1 maxK1 minK2 maxK2 Assembly_Type Number_of_Processors"
exit 1
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

ATYPE=$5
NUMProc=$6
ls *.F.fq.gz > namelist
sed -i'' -e 's/.F.fq.gz//g' namelist
NAMES=( `cat "namelist" `)

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
			pearRM -f $i.F.fq.gz -r $i.R.fq.gz -o $i -j $NUMProc -n $LENGTH &>kopt.log
			done
		cat namelist | parallel --no-notice -j $NUMProc "mawk '$AWK1' {}.assembled.fastq | mawk '$AWK2' | perl -e '$PERLT' > {}.uniq.seqs"
	fi
		
	
fi

#Create a data file with the number of unique sequences and the number of occurrences
if [ "$ATYPE" == "RPE" ]; then
	parallel --no-notice mawk -v x=$1 \''$1 >= x'\' ::: *.uniq.seqs | cut -f2 |  sort | uniq -c -w $FRL | sed -e 's/^[ \t]*//' | sed -e 's/\s/\t/g' | mawk -v x=$2 '$1 >= x' > uniq.k.$1.c.$2.seqs
else
	parallel --no-notice mawk -v x=$1 \''$1 >= x'\' ::: *.uniq.seqs | cut -f2 | perl -e 'while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}' | mawk -v x=$2 '$1 >= x' > uniq.k.$1.c.$2.seqs
fi
cut -f2 uniq.k.$1.c.$2.seqs > totaluniqseq
mawk '{c= c + 1; print ">dDocent_Contig_" c "\n" $1}' totaluniqseq > uniq.full.fasta
LENGTH=$(mawk '!/>/' uniq.full.fasta  | mawk '(NR==1||length<shortest){shortest=length} END {print shortest}')
LENGTH=$(($LENGTH * 3 / 4))
awk 'BEGIN {RS = ">" ; FS = "\n"} NR > 1 {print "@"$1"\n"$2"\n+""\n"gensub(/./, "I", "g", $2)}' uniq.full.fasta > uniq.fq
java -jar $TRIMMOMATIC SE -threads $NUMProc -phred33 uniq.fq uniq.fq1 ILLUMINACLIP:$ADAPTERS:2:30:10 MINLEN:$LENGTH &>/dev/null
mawk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' uniq.fq1 > uniq.fasta
mawk '!/>/' uniq.fasta > totaluniqseq
rm uniq.fq*

#If this is a PE assebmle
if [[ "$ATYPE" == "PE" || "$ATYPE" == "RPE" ]]; then
	#Reads are first clustered using only the Forward reads using CD-hit instead of rainbow
	sed -e 's/NNNNNNNNNN/\t/g' uniq.fasta | cut -f1 > uniq.F.fasta
	CDHIT=$(python -c "print max("$3" - 0.1,0.8)")
	cd-hit-est -i uniq.F.fasta -o xxx -c $CDHIT -T 0 -M 0 -g 1 -d 100 &>cdhit.log
	mawk '{if ($1 ~ /Cl/) clus = clus + 1; else  print $3 "\t" clus}' xxx.clstr | sed 's/[>dDocent_Contig_,...]//g' | sort -g -k1 > sort.contig.cluster.ids
	paste sort.contig.cluster.ids totaluniqseq > contig.cluster.totaluniqseq
	sort -k2,2 -g contig.cluster.totaluniqseq | sed -e 's/NNNNNNNNNN/\t/g' > rcluster
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
	sed -e 's/NNNNNNNNNN/\t/g' rainbow.fasta | cut -f1 | awk 'BEGIN {RS = ">" ; FS = "\n"} NR > 1 {print "@"$1"\n"$2"\n+""\n"gensub(/./, "I", "g", $2)}' > ref.F.fq
	sed -e 's/NNNNNNNNNN/\t/g' rainbow.fasta | cut -f2 | awk 'BEGIN {RS = ">" ; FS = "\n"} NR > 1 {print "@"$1"\n"$2"\n+""\n"gensub(/./, "I", "g", $2)}' > ref.R.fq

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
	paste other.F other.R | mawk '{if ($1 ~ />/) print $1; else print $0}' | sed 's/\t/NNNNNNNNNN/g' > other.FR

	cat other.FR overlap.fasta > totalover.fasta

	rm *.F *.R
fi
if [[ "$ATYPE" != "PE" && "$ATYPE" != "RPE" ]]; then
	cp uniq.fasta totalover.fasta
fi
cd-hit-est -i totalover.fasta -o reference.fasta.original -M 0 -T 0 -c $3 &>cdhit.log

sed -e 's/^C/NC/g' -e 's/^A/NA/g' -e 's/^G/NG/g' -e 's/^T/NT/g' -e 's/T$/TN/g' -e 's/A$/AN/g' -e 's/C$/CN/g' -e 's/G$/GN/g' reference.fasta.original > reference.fasta

SEQS=$(cat reference.fasta | wc -l)
SEQS=$(($SEQS / 2 ))
echo $SEQS

}

rm kopt.data &>/dev/null

for ((P = $1; P <= $2; P++))
	do
	for ((i = $3; i <= $4; i++))
	do
	X=$(($P + $i))
	if [ "$X" != "2" ]; then
	echo "K1 is $P" "K2 is $i" "c is 0.80"
	SEQS=$(Reference $P $i 0.8)
	echo $P $i 0.80 $SEQS >> kopt.data
		for j in {0.82,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98}
		do
		echo "K1 is $P" "K2 is $i" "c is $j"
		SEQS=$(Reference $P $i $j)
		echo $P $i $j $SEQS >> kopt.data
		done
	fi	
	done
done

cut -f4 -d " " kopt.data > plot.kopt.data
gnuplot << \EOF
set terminal dumb size 120, 30
set autoscale
unset label
set title "Histogram of number of reference contigs"
set ylabel "Number of Occurrences"
set xlabel "Number of reference contigs"
max = `sort -g plot.kopt.data | tail -1`
binwidth = max/250.0
bin(x,width)=width*floor(x/width) + binwidth/2.0
#set xtics 10
plot 'plot.kopt.data' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF

AF=$(mawk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' plot.kopt.data)
echo "Average contig number = $AF"
echo "The top three most common number of contigs"
echo -e "X\tContig number"
perl -e 'while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}' plot.kopt.data | sort -k1 -g -r | head -3
echo "The top three most common number of contigs (with values rounded)"
echo -e "X\tContig number"
while read NAME; do python -c "print round($NAME,-2)"; done < plot.kopt.data | perl -e 'while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}' | sort -g -r | head -3 | sed "s/^[ \t]*//"
