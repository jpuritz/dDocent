#!/bin/bash

#sed -i 's/\..*//g' $1

NAMES=( `cat $1 | cut -f1 | sort | uniq `)

for i in "${NAMES[@]}"
do
grep -A1 $i reference.fasta >> outliers.fasta 
done

