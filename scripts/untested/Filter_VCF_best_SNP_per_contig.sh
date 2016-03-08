#!/bin/bash

#Script to take the best SNP from every contig in a vcffile

if [[ -z "$1" ]]; then
echo "Usage is bash Filter_VCF_best_SNP_per_contig.sh vcffile"
exit 1
fi

NAME=$(echo $1 | sed -e 's/\.recode.*//g') 

cat $1 | mawk 'BEGIN{last_loc = 0} { 
		if ($1 ~/#/) print $0;
		else if ($1 ~!/#/ && last_loc == 0) {last_contig=$0; last_loc=$1; last_qual=$6}
		else if ($1 == last_loc && $6 > last_qual) {last_contig=$0; last_loc=$1; last_qual=$6}
		else if ($1 != last_loc) {print last_contig; last_contig=$0; last_loc=$1; last_qual=$6}
		} END{print last_contig}' > $NAME.filtered1SNPper.vcf


echo "Filtered VCF file is saved under name" $NAME.filtered1SNPper.vcf

