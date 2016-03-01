#!/usr/bin/env bash

if [ -z "$1" ]
then
echo "No file with barcodes and sample names specified."
echo "Correct usage: Rename_for_dDocent.sh barcodefile"
exit 1
else
NAMES=( `cut -f1  $1 `)
BARCODES=( `cut -f2 $1 `)
LEN=( `wc -l $1 `)
LEN=$(($LEN - 1))

echo ${NAMES[0]}
echo ${BARCODES[0]}

for ((i = 0; i <= $LEN; i++));
do
mv sample_${BARCODES[$i]}.1.fq.gz ${NAMES[$i]}.F.fq.gz
mv sample_${BARCODES[$i]}.2.fq.gz ${NAMES[$i]}.R.fq.gz
done
fi
