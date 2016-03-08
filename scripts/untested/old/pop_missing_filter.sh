#!/bin/bash

if [[ -z "$2" ]]; then
echo "Usage is pop_missing_filter vcffile popmap percent_missing_per_pop number_of_pops_for_cutoff name_for_output"
exit 1
fi

POPS=( `cut -f2 $2 | sort | uniq `)
rm badloci

for i in "${POPS[@]}"
do
grep -w $i $2 | cut -f1 > keep.$i
vcftools --vcf $1 --keep keep.$i --missing --out $i 
mawk '!/CHROM/' $i.lmiss | mawk -v x=$3 '$6 > x' | cut -f1,2 >> badloci
done

mawk '!/CH/' badloci | perl -e 'while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}' | mawk -v x=$4 '$1 > x' | cut -f2,3  > loci.to.remove

#sort badloci | uniq > loci.to.remove

vcftools --vcf $1 --exclude-positions loci.to.remove --recode --recode-INFO-all --out $5

