#!/bin/bash

echo "This script counts the number of potential genotyping errors due to low read depth, assuming that the probability of observing either allele in a heterozygote is 0.5"

R1H01S0=$(grep -wo 0/1:1:0 $1 | wc -l) 
R1H01S1=$(grep -wo 0/1:1:1 $1 | wc -l) 
R1H0S1=$(grep -wo 0/0:1:0 $1 | wc -l) 
R1H0S0=$(grep -wo 0/0:1:1 $1 | wc -l) 
R1H1S1=$(grep -wo 1/1:1:0 $1 | wc -l) 
R1H1S0=$(grep -wo 1/1:1:1 $1 | wc -l) 

R1GEN=$(python -c "print $R1H01S0+$R1H01S1+$R1H0S0+$R1H0S1+$R1H1S0+$R1H1S1")
R1ERR=$(python -c "print $R1GEN/2")

echo $R1ERR "maximum genotyping errors from genotypes from only 1 read"

R2H01S0=$(grep -wo 0/1:2:0 $1 | wc -l) 
R2H01S2=$(grep -wo 0/1:2:2 $1 | wc -l) 
R2H0S1=$(grep -wo 0/0:2:0 $1 | wc -l) 
R2H0S01=$(grep -wo 0/0:2:1 $1 | wc -l)
R2H0S0=$(grep -wo 0/0:2:2 $1 | wc -l)
R2H1S01=$(grep -wo 1/1:2:2 $1 | wc -l) 
R2H1S0=$(grep -wo 1/1:2:1 $1 | wc -l) 
R2H1S1=$(grep -wo 1/1:2:0 $1 | wc -l)

R2GEN=$(python -c "print $R2H01S0+$R2H01S1+$R2H0S0+$R2H0S1+$R2H0S01+$R2H1S0+$R2H1S1+$R2H1S01")
R2ERR=$(python -c "print $R2GEN/4")

echo $R2ERR "maximum genotyping errors from genotypes with only 2 reads"

R3H01S0=$(grep -wo 0/1:3:0 $1 | wc -l) 
R3H01S2=$(grep -wo 0/1:3:3 $1 | wc -l) 
R3H0S1=$(grep -wo 0/0:3:0 $1 | wc -l) 
R3H0S01=$(grep -wo 0/0:3:1 $1 | wc -l)
R3H0S0=$(grep -wo 0/0:3:2 $1 | wc -l)
R3H0S3=$(grep -wo 0/0:3:3 $1 | wc -l)
R3H1S01=$(grep -wo 1/1:3:2 $1 | wc -l) 
R3H1S0=$(grep -wo 1/1:3:1 $1 | wc -l) 
R3H1S1=$(grep -wo 1/1:3:0 $1 | wc -l)
R3H1S3=$(grep -wo 1/1:3:3 $1 | wc -l)

R3GEN=$(python -c "print $R3H01S0+$R3H01S1+$R3H0S0+$R3H0S1+$R3H0S01+$R3H1S0+$R3H1S1+$R3H1S01+$R3H0S3+$R3H1S3")
R3ERR=$(python -c "print $R3GEN/8")

echo $R3ERR "maximum genotyping errors from genotypes with only 3 reads"


IND=$(mawk '/#/' $1 | tail -1 | wc -w)
IND=$(($IND - 9))
LOCI=$(mawk '!/#/' $1 | wc -l)
GENO=$(( $IND * $LOCI ))

echo $IND "number of individuals and" $LOCI "equals" $GENO "total genotypes"

TOTERR=$(python -c "print $R1ERR+$R2ERR+$R3ERR")
ERRRATE=$(python -c "print $TOTERR/float($GENO)")

echo "Total potential error rate equals" $ERRRATE
