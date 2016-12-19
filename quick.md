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
  * This only needs to capture the total variablilty present in your data set
  * Include all `.fq.gz` files for each individual
7. Run [ReferenceOpt.sh](https://github.com/jpuritz/dDocent/blob/master/scripts/ReferenceOpt.sh)
  * Visualize data in `kopt.data`
    * Plot values for each k1,k2 combination across similarity thresholds
    * Pick a similarity threshold at the point of inflection on the curve
8. Run [RefMapOpt.sh](https://github.com/jpuritz/dDocent/blob/master/scripts/RefMapOpt.sh) using the similarity threshold picked from step 7.
  * Pick optimal k1,k2 cutoffs.  Ideally, you want to maximize properly paired mappings and coverage while minumizing mismatched reads.
9. Run dDocent on this subset with the correct assembly parameters, skipping mapping and snp calling.
10. Copy the `reference.fasta` file from this `RefOpt` directory to your main working directory.
11. Run dDocent on your full data set, skipping trimming and assembly.
12. Filter SNPs using what you've learned from [SNP Filtering Tutorial](/filtering), but also do your own exploring.
13. Repeat.
