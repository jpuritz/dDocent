#!/usr/bin/env bash
#dDocent needs an UTF-8 locale (gnuplot needs a point as the decimal
#separator), but exporting a locale that is not installed makes every
#program print setlocale warnings to stderr
if locale -a 2>/dev/null | grep -qiE '^en_us\.utf-?8$'; then
	export LC_ALL=en_US.UTF-8
fi

if [ -z "$2" ]; then
echo "Correct usage is sh remove.bad.hap.loci.sh file_with_bad_Loci vcf_file"
exit 1
fi

NAME=$(echo $2 | sed -e 's/\.recode.*//g') 

grep -vwf <(cut -f1 $1) $2 > $NAME.filtered.vcf
