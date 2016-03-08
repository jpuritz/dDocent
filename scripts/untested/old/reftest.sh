#!/bin/bash

rm kopt.data

for i in {3..30}
do
echo "K is $i"
SEQS=$(./reference.sh $i 2>/dev/null | mawk '/KSEQS/' | cut -f2)
echo $i $SEQS >> kopt.data
done


#Plot graph of above data
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Number of Unique Sequences with More than X Occurrences"
#set xlabel "Number of Occurrences"
set ylabel "Number of Unique Sequences
#set yr [0:100000]
plot 'kopt.data' with dots notitle
pause -1
EOF

