#!/bin/env bash

#This script will automatically remove sites in VCF files that do not have congruent genotypes across duplicate individuals
#It will automatically only consider genotypes that have at least 5 reads

#

if [[ -z "$2" ]]; then
echo "Usage is bash dup_sam_filter.sh VCF_file [File with duplicate sample names]"
echo "The list of names should have one line per pair of duplicate samples with tab separating the two names for the same individual"
exit 1
fi

echo "This script assumes that duplicate samples are named in the convention of PopA_001 and PopA_001a"

NAMES=( `cut -f1  $2 `)
NAM=( `cut -f2 $2 `)
LEN=( `wc -l $2 `)
LEN=$(($LEN - 1))


for ((i = 0; i <= $LEN; i++));
do
echo "${NAMES[$i]}" > keep.${NAMES[$i]}
echo "${NAM[$i]}" > keep.${NAM[$i]}

vcftools --vcf $1 --keep keep.${NAMES[$i]} --recode --recode-INFO-all --minDP 5 --out ${NAMES[$i]}
vcftools --vcf $1 --keep keep.${NAM[$i]} --recode --recode-INFO-all --minDP 5 --out ${NAM[$i]}

paste <(mawk '!/#/' ${NAMES[$i]}.recode.vcf | cut -f1,2,10 | cut -f1 -d ":") <(mawk '!/#/' ${NAM[$i]}.recode.vcf | cut -f1,2,10 | cut -f1 -d ":") | mawk '$3 != $6' | mawk '!/\./' | cut -f1,2 > bad.loci.${NAMES[$i]}.${NAM[$i]}

done


NAMES=( `cut -f1  $2 | sort | uniq `)
LEN=( `cut -f1 $2 | sort | uniq | wc -l `)
LEN=$(($LEN - 1))


cat bad.loci.${NAMES[0]}.* > total.bad.loci
rm ${NAMES[0]}.recode.vcf ${NAM[0]}.recode.vcf keep.${NAMES[0]} keep.${NAM[0]}

for ((i = 1; i <= $LEN; i++));
do
cat bad.loci.${NAMES[$i]}.* >> total.bad.loci 
rm ${NAMES[$i]}.recode.vcf ${NAM[$i]}.recode.vcf keep.${NAMES[$i]} keep.${NAM[$i]}

done

cat total.bad.loci | perl -e 'while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}' > mismatched.loci
rm total.bad.loci
