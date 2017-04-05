---
layout: page
title: Quick Start Guide
subtitle: Who reads manuals?
---

**Disclaimer** This is a work in progress


# Don't rush this thing.  Do it the right way!

1. Read the [User Guide](/UserGuide)!
2. Finish the [Assembly Tutorial](/assembly)!
3. Finish the [SNP Filtering Tutorial](/filtering)!
4. After properly demultiplexing and naming your data, run dDocent once with read trimming to test install.
5. Create a new folder called `RefOpt`
6. Place a subset of individuals of your total data set
  * Less than 5 individuals per region or locality
  * This only needs to capture the total variability present in your data set
  * Include all `.fq.gz` files for each individual
7. Run [ReferenceOpt.sh](https://github.com/jpuritz/dDocent/blob/master/scripts/ReferenceOpt.sh)
  * Visualize data in `kopt.data`
    * Plot values for each k1,k2 combination across similarity thresholds
    * Pick a similarity threshold at the point of inflection on the curve
  * Here's some R code to quickly plot this for you:

```R
library(ggplot2)

data.table <- read.table("kopt.data", header = FALSE, col.names= c("k1","k2","Similarity", "Contigs"))

data.table$K1K2 <- paste(data.table$k1, data.table$k2, sep=",")

df=data.frame(data.table)
df$Kcombo <- as.factor(df$K1K2)

p <- ggplot(df, aes(x=Similarity, y=Contigs, group=K1K2)) + scale_x_continuous(breaks=seq(0.8,0.98,0.01)) + geom_line(aes(colour = K1K2))
p
```
       
8. Run [RefMapOpt.sh](https://github.com/jpuritz/dDocent/blob/master/scripts/RefMapOpt.sh) using the similarity threshold picked from step 7. 
  * **Note- You will need to have the trimmed reads files `*.R1.fq.gz` and `*.R2.fq.gz` included to run this script**
  * Pick optimal k1,k2 cutoffs.  Ideally, you want to maximize properly paired mappings and coverage while minimizing mismatched reads.
9. Run dDocent on this subset with the correct assembly parameters, skipping mapping and snp calling.
10. Copy the `reference.fasta` file from this `RefOpt` directory to your main working directory.
11. Run dDocent on your full data set, skipping trimming and assembly.
12. Filter SNPs using what you've learned from [SNP Filtering Tutorial](/filtering), but also do your own exploring.
13. Repeat.
