#!/bin/bash


NAME=$(echo $2 | sed -e 's/\.recode.*//g') 

grep -vwf <(cut -f1 $1) $2 > $NAME.filtered.vcf

