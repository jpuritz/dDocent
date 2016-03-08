#!/bin/bash

vcftools --vcf SNP-VQSR-PASS.vcf --minQ 20 --geno 0.5 --mac 2 --recode --recode-INFO-all --out TRSg
vcftools --vcf TRSg.recode.vcf --minDP 4 --recode --recode-INFO-all --out TRSgdp4 &
vcftools --vcf TRSg.recode.vcf --minDP 3 --recode --recode-INFO-all --out TRSgdp3 &
vcftools --vcf TRSg.recode.vcf --minDP 5	--recode --recode-INFO-all --out TRSgdp5 &
vcftools --vcf TRSg.recode.vcf --minDP 10 --recode --recode-INFO-all --out TRSgdp10 &


