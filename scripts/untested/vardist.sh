#!/bin/bash

vcfdistance <$1 | grep -oh 'BasesToClosestVariant=[0-9]*' | sed -e 's/BasesToClosestVariant=//g' > vardist

gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale
#set xrange [10:150]
unset label
set title "Histogram of distance to closest variant"
set ylabel "Number of Occurrences"
set xlabel "Mean Depth"
#set yr [0:100000]
binwidth=1
bin(x,width)=width*floor(x/width) + binwidth/2.0
set xtics 5
plot 'vardist' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF


