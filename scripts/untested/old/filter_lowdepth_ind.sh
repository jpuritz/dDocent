#!/bin/bash

vcftools --vcf $1 --depth --out $2

CUTOFF=$(mawk '!/IN/' $2.idepth | cut -f3 | sort -rn | perl -e '$d=.85;@l=<>;print $l[int($d*$#l)]')
mawk -v x=$CUTOFF '$3 < x' $2.idepth | cut -f1 > lowDP.indv

vcftools --vcf $1 --remove lowDP.indv --recode --recode-INFO-all --out $2
