#!/bin/bash

#sed -i 's/\..*//g' $1

NAMES=( `mawk '!/#/' $1 | cut -f1 | sort | uniq `)

rm outliers.fasta

for i in "${NAMES[@]}"
do
grep -wA1 $i reference.fasta >> outliers.fasta 
done

