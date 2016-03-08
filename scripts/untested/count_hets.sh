#!/bin/bash
vcftools --vcf $1 --max-alleles 2 --recode --out out --site-depth &>/dev/null
HET=( `mawk '!/#/ {print $1,gsub(/0[\/\|]1/,"")}' out.recode.vcf | cut -f2 -d " "` )
HETT=( `mawk '!/#/ {print $1,gsub(/1[\/\|]0/,"")}' out.recode.vcf | cut -f2 -d " "` )
DEPTH=( `mawk '!/SUM/' "out.ldepth" | cut -f3`)

LEN=${#HET[@]}
LEN=$((LEN - 1))
rm hapcounts 2>/dev/null

for ((i = 0; i <= $LEN; i++));
do
HETTT=$((${HET[$i]} + ${HETT[$i]} ))
echo $HETTT ${DEPTH[$i]} >> hapcounts
done

gnuplot << \EOF
set terminal dumb size 120, 30
set autoscale
f(x) = a * x + b
fit f(x) "hapcounts" using 2:1 via b, a
unset label
set title "Number of Unique Sequences with More than X Occurrences"
set xlabel "Number of Occurrences"
set ylabel "Number of Unique Sequences
plot 'hapcounts' using 2:1 with dots, f(x) title "Model Fit" with points
pause -1
EOF
