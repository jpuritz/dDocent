#!/bin/bash

ls *-RG.bam >bamlist.list; bamtools merge -list bamlist.list > full.bam
samtools index full.bam &>/dev/null
samtools idxstats full.bam 2>/dev/null | cut -f3 > coverage
rm full.bam*

gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale
#set xrange [10:150] 
unset label
set title "Histogram of coverage per contig"
set ylabel "Number of Occurrences"
set xlabel "Coverage"
#set yr [0:100000]
binwidth=10
bin(x,width)=width*floor(x/width) + binwidth/2.0
set xtics 100
plot 'coverage' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF

AVE=$(mawk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' coverage)
NZ=$(mawk '$1 == 0' coverage | wc -l)
NZ=$(($NZ + 0))

echo "The average coverage per contig is" $AVE
echo "There are" $NZ "contigs with zero coverage"
