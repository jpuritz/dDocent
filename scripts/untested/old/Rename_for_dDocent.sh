#!/bin/bash

if [ -z "$1"]
then
echo "No file with barcodes and sample names specified."
echo "Correct usage: Rename_for_dDocent.sh barcodefile"
exit 1
else
NAMES=( `cut -f1  $1 `)
BARCODES=( `cut -f2 $1 `)
LEN=( `wc -l $1 `)
LEN=$(($LEN - 0))

echo ${NAMES[0]}
echo ${BARCODES[0]}

for ((i = 0; i <= $LEN; i++));
do
mv sample_${BARCODES[$i]}.1.fq ${NAMES[$i]}.F.fq
mv sample_${BARCODES[$i]}.2.fq ${NAMES[$i]}.R.fq
done
fi

