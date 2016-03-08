#!/bin/bash

vcftools --vcf $1 --missing --out $2

CUTOFF=$(mawk '!/IN/' $2.imiss | cut -f5 | sort -rn | perl -e '$d=.14;@l=<>;print $l[int($d*$#l)]')
#echo $CUTOFF

mawk '!/IN/' $2.imiss | cut -f5 > totalmissing

gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Histogram of % missing data per individual"
set ylabel "Number of Occurrences"
set xlabel "% of missing data"
#set yr [0:100000]
binwidth=0.01
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'totalmissing' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF

echo "The 85% cutoff would be" $CUTOFF
echo "Would you like to set a different cutoff, yes or no"

read NEWCUTOFF

if [ "$NEWCUTOFF" != "yes" ]; then

mawk -v x=$CUTOFF '$5 > x' $2.imiss | cut -f1 > lowDP.indv

vcftools --vcf $1 --remove lowDP.indv --recode --recode-INFO-all --out $2

else

echo "Please enter new cutoff"

read CUTOFF2

mawk -v x=$CUTOFF2 '$5 > x' $2.imiss | cut -f1 > lowDP.indv

vcftools --vcf $1 --remove lowDP.indv --recode --recode-INFO-all --out $2
fi
