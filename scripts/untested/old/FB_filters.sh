#!/bin/bash

#!/bin/bash

if [[ -z "$2" ]]; then
echo "Usage is sh FB_filters.sh VCF_file Output_prefix"
exit 1
fi

IND=$(mawk '/#/' $1 | tail -1 | wc -w)
IND=$(($IND - 9))
vcffilter -f "AB > 0.28" $1 > $2
vcffilter -f "SAR > 50 & SAF < 50 | SAF > 50 & SAR < 50" -t PASS -F OL $2 > $2.fil1.vcf
vcffilter -f "PAIRED < 0.1 & PAIREDR > 0.1 | PAIRED > 0.1 & PAIREDR < 0.1" -t NP -F PASS -A $2.fil1.vcf > $2.fil2.vcf
vcffilter -f "PAIRED > 0.75" -t PASS -F NP2 -A $2.fil2.vcf > $2.fil3.vcf
vcffilter -f "SRR > 50 & SRF < 50 | SRF > 50 & SRR < 50" -t PASS -F OL2 -A $2.fil3.vcf > $2.fil4.vcf
vcftools --vcf $2.fil4.vcf --keep-filtered PASS,PASS,PASS,PASS --keep-filtered PASS,PASS,NP2,PASS --site-depth --out $1 
cut -f3 $1.ldepth > $1.site.depth
DP=$(mawk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' $1.site.depth)
#DP=$(mawk '{print $1}' $1.site.depth | sort -rn | perl -e '$d=.025;@l=<>;print $l[int($d*$#l)]')
SD=$(mawk '{delta = $1 - avg; avg += delta / NR; mean2 += delta * ($1 - avg); } END { print sqrt(mean2 / NR); }' $1.site.depth)
#SD=$(awk '{sum+=$1; sumsq+=$1*$1} END {print sqrt(sumsq/NR - (sum/NR)**2)}' $1.site.depth)
DP=$(python -c "print ($DP+ $SD) / $IND")
echo $DP
vcftools --vcf $2.fil4.vcf --keep-filtered PASS,PASS,PASS,PASS --keep-filtered PASS,PASS,NP2,PASS --recode-INFO-all --out $2.FIL --max-meanDP $DP --recode 
