#!/bin/bash

if [[ -z "$1" ]]; then
echo "Usage: LositantoOutlierVCF.sh Vcffile LositanFile FDR Prefix_for_output"
exit 1
fi

PVALUE=$(echo "1-$3/2" | bc -l)

mawk '!/#/' $1 | cut -f1,2 > totalloci
mawk '!/Het/' $2 | cut -f2-5 > LS.noheader
paste totalloci LS.noheader | awk -v x=$PVALUE '$5 >= x' | cut -f1,2 > Outlier.list
vcftools --vcf $1 --positions Outlier.list --recode --recode-INFO-all --out $4.outlieronly
vcftools --vcf $1 --exclude-positions Outlier.list --recode --recode-INFO-all --out $4.neutralonly

rm totalloci
rm LS.noheader
mv Outlier.list $4.Outlier.list
