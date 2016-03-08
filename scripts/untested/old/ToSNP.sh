#!/bin/bash

NAME=$(echo $1 | sed -e 's/\.recode.*//g')

vcfallelicprimitives --keep-info --keep-geno $1 > $NAME.prim.vcf
vcftools --vcf $NAME.prim.vcf --remove-indels --recode --recode-INFO-all --out SNP.$NAME
