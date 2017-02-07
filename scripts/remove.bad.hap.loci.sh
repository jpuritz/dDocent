#!/usr/bin/env bash
export LC_ALL=en_US.UTF-8

if [ -z "$2" ]; then
echo "Correct usage is sh remove.bad.hap.loci.sh file_with_bad_Loci vcf_file"
exit 1
fi

NAME=$(echo $2 | sed -e 's/\.recode.*//g') 

grep -vwf <(cut -f1 $1) $2 > $NAME.filtered.vcf
