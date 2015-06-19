#!/bin/bash

##THIS NEEDS AN UPDATE

if [[ -z "$2" ]]; then
echo "Usage is sh ReferenceOpt.sh K_value c_value"
exit 1
fi

if [ -e "uniq.seqs" ]; then
i=0
  else
    echo "This script needs to be run after an initial dDocent assembly."
    exit 1
    fi


if [ -f "uniq.seqs.gz" ]; then
	if [ uniq.seqs.gz -nt uniq.seqs ]; then
	gunzip uniq.seqs.gz 2>/dev/null
	fi
fi

mawk -v x=$1 '$1 >= x' uniq.seqs | cut -f 2 > totaluniqseq
rm uniq.fasta &>/dev/null
cat totaluniqseq | while read line
	do
	echo ">Contig"$i >>uniq.fasta
	echo $line >> uniq.fasta
	i=$(($i + 1))
	done

sed -e 's/NNNNNNNNNN/\t/g' uniq.fasta | cut -f1 > uniq.F.fasta
sed -e 's/NNNNNNNNNN/\t/g' uniq.fasta | cut -f2 > uniq.R.fasta

seqtk seq -r uniq.R.fasta > uniq.RC.fasta
rm uniq.R.fasta

#Now use rainbow to cluster and assemble reads into longer contigs
rainbow cluster -m 6 -1 uniq.F.fasta -2 uniq.RC.fasta > rcluster 2> rainbow.log
rainbow div -i rcluster -o rbdiv.out -f 0.01
rainbow merge -o rbasm.out -a -i rbdiv.out -r 2
select_best_rbcontig_plus_read1.pl rbasm.out rbdiv.out >rainbow.fasta

#cd-hit to cluster reads based on sequence similarity
cd-hit-est -i rainbow.fasta -o referenceRC.fasta -M 0 -T 0 -c $2 &>cdhit.log

seqtk seq -r referenceRC.fasta > reference.fasta.original
rm referenceRC.fasta

sed -e 's/^C/NC/g' -e 's/^A/NA/g' -e 's/^G/NG/g' -e 's/^T/NT/g' -e 's/T$/TN/g' -e 's/A$/AN/g' -e 's/C$/CN/g' -e 's/G$/GN/g' reference.fasta.original > reference.fasta

SEQS=$(cat reference.fasta | wc -l)
SEQS=$(($SEQS / 2 ))
echo $SEQS

samtools faidx reference.fasta &>/dev/null
bwa index reference.fasta &>/dev/null
