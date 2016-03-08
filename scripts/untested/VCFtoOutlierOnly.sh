#!/bin/bash

if [[ -z "$1" ]]; then
echo "Usage: VCFtoOutlierOnly.sh Vcffile BayescanOutput FDR Prefix_for_output"
exit 1
fi

mawk '!/#/' $1 | cut -f1,2 > totalloci
mawk '!/qval/' $2 | cut -f2-6 > BS.noheader
paste totalloci BS.noheader | mawk -v x=$3 '$6 <= x' | cut -f1,2 > Outlier.list
vcftools --vcf $1 --positions Outlier.list --recode --recode-INFO-all --out $4.outlieronly
vcftools --vcf $1 --exclude-positions Outlier.list --recode --recode-INFO-all --out $4.neutralonly
