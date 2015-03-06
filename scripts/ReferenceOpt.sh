#!/bin/bash

if [[ -z "$2" ]]; then
echo "Usage is sh ReferenceOpt.sh lowK highK"
exit 1
fi

Reference(){


if [ -f "uniq.seqs.gz" ]; then
	if [ uniq.seqs.gz -nt uniq.seqs ]; then
	gunzip uniq.seqs.gz 2>/dev/null
	fi
fi

mawk -v x=$1 '$1 >= x' uniq.seqs | cut -f 2 > totaluniqseq
rm uniq.fasta &>/dev/null
cat totaluniqseq | while read line
	do
	echo ">Contig"$i >>uniq.fasta
	echo $line >> uniq.fasta
	i=$(($i + 1))
	done

sed -e 's/NNNNNNNNNN/\t/g' uniq.fasta | cut -f1 > uniq.F.fasta
sed -e 's/NNNNNNNNNN/\t/g' uniq.fasta | cut -f2 > uniq.R.fasta

seqtk seq -r uniq.R.fasta > uniq.RC.fasta
rm uniq.R.fasta

#Now use rainbow to cluster and assemble reads into longer contigs
rainbow cluster -m 6 -1 uniq.F.fasta -2 uniq.RC.fasta > rcluster 2> rainbow.log
rainbow div -i rcluster -o rbdiv.out -f 0.01
rainbow merge -o rbasm.out -a -i rbdiv.out -r 2
select_best_rbcontig_plus_read1.pl rbasm.out rbdiv.out >rainbow.fasta

#cd-hit to cluster reads based on sequence similarity
cd-hit-est -i rainbow.fasta -o reference.fasta -M 0 -T 0 -c $2 &>cdhit.log

SEQS=$(cat reference.fasta | wc -l)
SEQS=$(($SEQS / 2 ))
echo $SEQS
}

rm kopt.data &>/dev/null

for ((i = $1; i <= $2; i++))
do
echo "K is $i" "c is 0.80"
SEQS=$(Reference $i 0.80)
echo $i 0.80 $SEQS >> kopt.data
	for j in {0.82,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98}
	do
	echo "K is $i" "c is $j"
	cd-hit-est -i rainbow.fasta -o reference.fasta -M 0 -T 0 -c $j &>cdhit.log
	SEQS=$(cat reference.fasta | wc -l)
	SEQS=$(($SEQS / 2 ))
	echo $i $j $SEQS >> kopt.data
	done
done

cut -f3 -d " " kopt.data > plot.kopt.data
gnuplot << \EOF
set terminal dumb size 120, 30
set autoscale
unset label
set title "Histogram of number of reference contigs"
set ylabel "Number of Occurrences"
set xlabel "Number of reference contigs"
max = `sort -g plot.kopt.data | tail -1`
binwidth = max/150.0
bin(x,width)=width*floor(x/width) + binwidth/2.0
#set xtics 10
plot 'plot.kopt.data' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF

AF=$(mawk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' plot.kopt.data)
echo "Average contig number = $AF"
echo "The top three most common number of contigs"
echo -e "X\tContig number"
perl -e 'while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}' plot.kopt.data | sort -k1 -g -r | head -3
echo "The top three most common number of contigs (with values rounded)"
echo -e "X\tContig number"
while read NAME; do python -c "print round($NAME,-2)"; done < plot.kopt.data | perl -e 'while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}' | sort -g -r | head -3 | sed "s/^[ \t]*//"
