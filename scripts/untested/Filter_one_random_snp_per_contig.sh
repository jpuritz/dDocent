#!/usr/bin/env bash

#Script to take random SNP from every contig in a vcffile

if [[ -z "$1" ]]; then
echo "Usage is bash Filter_one_random_snp_per_contig.sh vcffile"
exit 1
fi

#Calculate number of SNPs
Loci=(`mawk '!/#/' $1 | wc -l `)

#Generate list of random numbers
seq 1 500000 | shuf | head -$Loci > nq

#create temporary file that has a random number assigned to each SNP in first column
cat <(mawk '/^#/' $1) <(paste <(mawk '!/#/' $1 | cut -f1-5) nq <(mawk '!/#/' $1 | cut -f7- ) )> temp

#Get name of VCF file
NAME=$(echo $1 | sed -e 's/\.recode.*//g' | sed -e 's/.vcf//g' ) 

#Use awk (mawk) to parse file and select one snp per contig (one with largest random number)
cat temp | mawk 'BEGIN{last_loc = 0} { 
		if ($1 ~/#/) print $0;
		else if ($1 ~!/#/ && last_loc == 0) {last_contig=$0; last_loc=$1; last_qual=$6}
		else if ($1 == last_loc && $6 > last_qual) {last_contig=$0; last_loc=$1; last_qual=$6}
		else if ($1 != last_loc) {print last_contig; last_contig=$0; last_loc=$1; last_qual=$6}
		} END{print last_contig}' | mawk 'NF > 0' > $NAME.filtered1SNPper.vcf

#Remove temp file
rm temp

#Announce triumphant completion
echo "Filtered VCF file is saved under name" $NAME.filtered1SNPper.vcf
