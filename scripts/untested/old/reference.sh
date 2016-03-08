#!/bin/bash

NAMES=( `cat "namelist" `)


mawk -v x=$1 '$1 >= x' uniq.seqs | cut -f 2 > totaluniqseq


#Convert reads to fasta
uniq2fasta totaluniqseq > uniq.fasta

#Perl function to split contigs by length

sed -e 's/NNNNNNNNNN/\t/g' uniq.fasta | cut -f1 > uniq.F.fasta
sed -e 's/NNNNNNNNNN/\t/g' uniq.fasta | cut -f2 > uniq.R.fasta

#sed -i'' -e 's/_.* (.*)/_1/g' uniq.F.fasta
#sed -i'' -e 's/_.*_.* (.*)/_2/g' uniq.R.fasta 

seqtk seq -r uniq.R.fasta > uniq.RC.fasta
rm uniq.R.fasta

#Now use rainbow to cluster and assemble reads into longer contigs
rainbow cluster -m 6 -1 uniq.F.fasta -2 uniq.RC.fasta > rcluster
rainbow div -i rcluster -o rbdiv.out -f 0.01
rainbow merge -o rbasm.out -a -i rbdiv.out -r 2
select_best_rbcontig_plus_read1.pl rbasm.out rbdiv.out >rainbow.fasta

#cd-hit to cluster reads based on sequence similarity
cd-hit-est -i rainbow.fasta -o referenceRC.fasta -mask N -M 0 -T 0 -c 0.9 &>cdhit.log

seqtk seq -r referenceRC.fasta > reference.fasta.original

sed -e 's/^C/NC/g' -e 's/^A/NA/g' -e 's/^G/NG/g' -e 's/^T/NT/g' -e 's/T$/TN/g' -e 's/A$/AN/g' -e 's/C$/CN/g' -e 's/G$/GN/g' reference.fasta.original > reference.fasta

samtools faidx reference.fasta
bwa index reference.fasta

SEQS=$(cat reference.fasta | wc -l)
SEQS=$(($SEQS / 2 ))

echo -e "KSEQS\t" $SEQS
