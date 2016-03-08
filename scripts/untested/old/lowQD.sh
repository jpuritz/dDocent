#!/bin/bash


cut -f8 $1 | grep -oe "DP=[0-9]*" | sed -s 's/DP=//g' > $1.DEPTH
mawk '!/#/' $1 | cut -f1,2,6 > $1.loci.qual
DEPTH=$(mawk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' $1.DEPTH)
QUAL=$(python -c "print $DEPTH * 2")
paste $1.loci.qual $1.DEPTH | mawk -v x=$DEPTH '$4 > x'| mawk '$3 < 1.0 * $4' > $1.lowQDloci

#DEPTH=$(awk '{print $1}' $1.DEPTH | sort -rn | perl -e '$d=.25;@l=<>;print $l[int($d*$#l)]')
DEPTH=$(mawk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' $1.DEPTH)
SD=$(mawk '{delta = $1 - avg; avg += delta / NR; mean2 += delta * ($1 - avg); } END { print sqrt(mean2 / NR); }' $1.DEPTH)
DEPTH=$(python -c "print int("$DEPTH") + int("$SD")")
paste $1.loci.qual $1.DEPTH | mawk -v x=$DEPTH '$4 > x'| mawk '$3 < 2 * $4' >> $1.lowQDloci

vcftools --vcf $1 --exclude-positions $1.lowQDloci --recode --recode-INFO-all --out $1.LowQDFIL
echo $DEPTH
echo $QUAL
#vcffilter -f "DP > $DEPTH & QUAL < $QUAL" $1.inter.recode.vcf -t lowQD -F PASS > $1.inter2
#vcftools --vcf $1.inter2 --remove-filtered lowQD --recode --recode-INFO-all --out $1.LQDFIL
