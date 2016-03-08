#!/bin/bash

paste <(mawk '/Contig.*:/' $1 ) <(mawk '/Over/' $1) | mawk '$4 > 0.95' | cut -f1 | sed 's/://g' > Hap.loc.belowMAF05

LOC=$(mawk '/Contig/' $2 | wc -l)   ####This gives you the number of loci, use it in the next line of code

paste <(seq $LOC ) <(mawk '/Contig/' $2) > BS.loci.numbered

grep -wf Hap.loc.belowMAF05 <(sed 's/://g' BS.loci.numbered) | cut -f1 > bs.loci.lowmaf  

