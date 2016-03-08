#!/bin/bash

#Script to randomly sample 18 individuals from populations and calculate pi

mawk1='{ sum += $5; n++ } END { if (n > 0) print sum / n; }'

rm sample.pi sample.pi.pop subsample.pi

for i in {1..100}
do
mawk '/AR_/' kept.popmap | shuf | head -15 | cut -f1 > AR1.keep
mawk '/ARB_/' kept.popmap | shuf | head -15 | cut -f1 > AR2.keep
mawk '/MB_/' kept.popmap | shuf | head -15 | cut -f1 > NR1.keep
mawk '/MBB_/' kept.popmap | shuf | head -15 | cut -f1 > NR2.keep
mawk '$2 =="EL"' kept.popmap | shuf | head -15 | cut -f1 > NR4.keep
mawk '$2 =="ELA"' kept.popmap | shuf | head -15 | cut -f1 > NR3.keep
mawk '$2 =="JC"' kept.popmap | shuf | head -15 | cut -f1 > AP1.keep
mawk '$2 =="PC"' kept.popmap | shuf | head -15 | cut -f1 > AP2.keep
ls *.keep | parallel --no-notice vcftools --vcf SNP.TRSdp5MIp25g9HWEHFv2.neutralmaf025.recode.vcf --keep {} --out {.} --site-pi &> /dev/null
ls *.sites.pi | parallel --no-notice mv {} {.}
ls *.sites | parallel --no-notice "mawk '$mawk1' {} >> sample.pi && echo {.} >> sample.pi.pop"
done

paste sample.pi sample.pi.pop > subsample.pi
