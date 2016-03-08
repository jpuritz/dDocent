#!/bin/bash

NAME=$(echo $1 | sed -e 's/\.recode.*//g') 

vcfallelicprimitives --keep-info --keep-geno $1 > $NAME.prim.vcf
vcftools --vcf $NAME.prim.vcf --remove-indels --recode --recode-INFO-all --out SNP.$NAME >/dev/null

bgzip SNP.$NAME.recode.vcf
tabix -p vcf SNP.$NAME.recode.vcf.gz
vcf-annotate --filter c=3,25 SNP.$NAME.recode.vcf.gz > SNP.$NAME.vcf
vcftools --vcf SNP.$NAME.vcf --keep-filtered SnpCluster --recode --recode-INFO-all --out SNP.$NAME.1
mawk '!/#/' SNP.$NAME.1.recode.vcf | cut -f1,2 > SNP.$NAME.LOC
HET=( `mawk '!/#/ {print $1,gsub(/0[\/\|]1/,"")}' SNP.$NAME.1.recode.vcf  | cut -f2 -d " "` )
HETT=( `mawk '!/#/ {print $1,gsub(/1[\/\|]0/,"")}' SNP.$NAME.1.recode.vcf  | cut -f2 -d " "` )

LEN=${#HET[@]}
LEN=$((LEN - 1))
rm SNP.$NAME.hets 2>/dev/null

for ((i = 0; i <= $LEN; i++));
do
HETTT=$((${HET[$i]} + ${HETT[$i]} ))
echo $HETTT >> SNP.$NAME.hets
done

paste SNP.$NAME.LOC SNP.$NAME.hets > SNP.$NAME.hhets
CUT=$(cut -f3 SNP.$NAME.hhets | sort -g | perl -e '$d=.9;@l=<>;print $l[int($d*$#l)]')
mawk -v x=$CUT '$3 > x' SNP.$NAME.hhets | cut -f1,2 > snpclusters.to.filter

vcftools --gzvcf SNP.$NAME.recode.vcf.gz --exclude-positions snpclusters.to.filter --recode --recode-INFO-all --out SNP.$NAME.SCfiltered
if [[ -n "$2" ]]; then
vcftools --vcf SNP.$NAME.SCfiltered.recode.vcf --exclude-positions $2 --recode --recode-INFO-all --out SNP.$NAME.SCfilteredF
fi

#gnuplot << \EOF 
#set terminal dumb size 120, 30
#set autoscale
#set xrange [1:150]
#unset label
#set title "Histogram of number of heterozygotes"
#set ylabel "Number of Occurrences"
#set xlabel "Mean Depth"
#set yr [0:100000]
#binwidth=1
#bin(x,width)=width*floor(x/width) + binwidth/2.0
#set xtics 5
#plot 'SNP.$NAME.hhets' using (bin($1,binwidth)):(1.0) smooth freq with boxes
#pause -1
#EOF
