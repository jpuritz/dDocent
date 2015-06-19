#!/bin/bash

NAMES=( `cat $1 | cut -f1 | sort | uniq `)

NAME=$(echo $2 | sed -e 's/\.recode.*//g') 
LEN=${#NAMES[@]}

mawk -v x="${NAMES[0]}" '$0 !~ x' $NAME.recode.vcf > $NAME.0.t.vcf

for ((i = 1; i < $LEN; i++));
do
j=$(($i - 1))
mawk -v x=${NAMES[$i]} '$0 !~ x' $NAME.$j.t.vcf > $NAME.$i.t.vcf
done

LAST=$(($LEN - 1))

mv $NAME.$LAST.t.vcf $NAME.filtered.vcf

rm $NAME.*.t.vcf
