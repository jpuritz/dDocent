#!/usr/bin/env bash

if [[ -z "$3" ]]; then
echo "Usuage is multi.maf.sh [vcffile] [maf] [prefix for outfile]"
exit 1
fi


paste <(cut -f1,2 $1 | mawk '!/#/' ) <(cut -f8 $1 | grep -oh "AF=.*;AN" | sed -e 's/AF=//g' -e 's/;AN//g' -e 's/,/\t/g' | mawk '{ for(i=1; i<=NF;i++) j+=$i; print j; j=0 }') > all.maf.frq


mawk -v x=$2 '$3 > x && $3 < 1-x' all.maf.frq | cut -f1,2 > maf.loci.to.keep

PREFIX=$3

vcftools --vcf $1 --positions maf.loci.to.keep --recode --recode-INFO-all --out $PREFIX
