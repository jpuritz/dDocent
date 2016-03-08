#!/bin/bash

for i in {3..40}
do
echo "K is $i" >>test
sh ./reference.sh $i
bwa mem reference.fasta JC_1161.R1.fq JC_1161.R2.fq -L 100,5 -t 32 -a -M -T 10 -A 1 -B 3 -O 5 -R "@RG\tID:test\tSM:test\tPL:Illumina" 2> /dev/null | mawk '!/\t[2-9].S.*/' | mawk '!/[2-9].S\t/' | samtools view -@32 -q 1 -SbT reference.fasta - 2>> test | samtools flagstat - 2> /dev/null >> test
done


