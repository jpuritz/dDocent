#!/usr/bin/env bash

echo "This script counts the number of potential genotyping errors due to low read depth"
echo "It report a low range, based on a 50% binomial probability of observing the second allele in a heterozygote and a high range based on a 25% probability."

R1H01S0=$(grep -oh '0[/|]1:1:' $1 | wc -l) 
R1H01S1=$(grep -oh '1[/|]0:1:' $1 | wc -l) 
R1H0S1=$(grep -oh '0[/|]0:1:' $1 | wc -l) 
R1H1S1=$(grep -oh '1[/|]1:1:' $1 | wc -l) 

R1GEN=$(python -c "print $R1H01S0+$R1H01S1+$R1H0S1+$R1H1S1")
R1ERR1=$(python -c "print $R1GEN/2")
R1ERR2=$(python -c "print $R1GEN*0.75")

echo "Potential genotyping errors from genotypes from only 1 read range from $R1ERR1 to $R1ERR2" 

R2H01S0=$(grep -oh '0[/|]1:2:[02]' $1 | wc -l) 
R2H01S2=$(grep -oh '1[/|]0:2:[02]' $1 | wc -l) 
R2H0S1=$(grep -oh '0[/|]0:2:' $1 | wc -l) 
R2H1S1=$(grep -oh '1[/|]1:2:' $1 | wc -l)

R2GEN=$(python -c "print $R2H01S0+$R2H01S2+$R2H1S1+$R2H0S1")
R2ERR1=$(python -c "print $R2GEN/4")
R2ERR2=$(python -c "print $R2GEN* 0.5625")

echo "Potential genotyping errors from genotypes from only 2 reads range from $R2ERR1 to $R2ERR2" 

R3H01S0=$(grep -oh '0[/|]1:3:[03]' $1 | wc -l) 
R3H01S2=$(grep -oh '1[/|]0:3:[03]' $1 | wc -l) 
R3H0S1=$(grep -oh '0[/|]0:3:' $1 | wc -l)
R3H1S3=$(grep -oh '1[/|]1:3:' $1 | wc -l)

R3GEN=$(python -c "print $R3H01S0+$R3H01S2+$R3H0S1+$R3H1S3")
R3ERR1=$(python -c "print $R3GEN/8")
R3ERR2=$(python -c "print $R3GEN*0.42")

echo "Potential genotyping errors from genotypes from only 3 reads range from $R3ERR1 to $R3ERR2"

R4H0=$(grep -oh '0[/|]0:4:' $1 | wc -l)
R4H1=$(grep -oh '1[/|]1:4:' $1 | wc -l)

R4GEN=$(python -c "print $R4H0+$R4H1")
R4ERR1=$(python -c "print $R4GEN/16")
R4ERR2=$(python -c "print $R4GEN*0.316")

echo "Potential genotyping errors from genotypes from only 4 reads range from $R4ERR1 to $R4ERR2"

R5H0=$(grep -oh '0[/|]0:5:' $1 | wc -l)
R5H1=$(grep -oh '1[/|]1:5:' $1 | wc -l)

R5GEN=$(python -c "print $R5H0+$R5H1")
R5ERR1=$(python -c "print $R5GEN/32")
R5ERR2=$(python -c "print int(	$R5GEN*0.237)")

echo "Potential genotyping errors from genotypes from only 5 reads range from $R5ERR1 to $R5ERR2"

IND=$(mawk '/#/' $1 | tail -1 | wc -w)
IND=$(($IND - 9))
LOCI=$(mawk '!/#/' $1 | wc -l)
MISSING=$(grep -Fwo ./.:. $1 | wc -l)
GENO=$(( $IND * $LOCI  ))

echo $IND "number of individuals and" $LOCI "equals" $GENO "total genotypes"
GENO=$(( $GENO - $MISSING))
echo Total genotypes not counting missing data $GENO

TOTERR1=$(python -c "print $R1ERR1+$R2ERR1+$R3ERR1+$R4ERR1+$R5ERR1")
TOTERR2=$(python -c "print $R1ERR2+$R2ERR2+$R3ERR2+$R4ERR2+$R5ERR2")
ERRRATEL=$(python -c "print $TOTERR1/float($GENO)")
ERRRATEH=$(python -c "print $TOTERR2/float($GENO)")

echo "Total potential error rate is between $ERRRATEL and $ERRRATEH"

ALL=$(($R1GEN+$R2GEN+$R3GEN+$R4GEN+$R5GEN))
ERRALL=$(python -c "print $ALL/float($GENO)")

echo "SCORCHED EARTH SCENARIO"
echo "WHAT IF ALL LOW DEPTH HOMOZYGOTE GENOTYPES ARE ERRORS?????"
echo "The total SCORCHED EARTH error rate is $ERRALL."
