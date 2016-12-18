---
layout: post
title: New Version of Pipeline **Important note for all previous versions**
tags: [From Old Website]
---

A never version of the pipeline is now up at the repository.  This version helps with OS X compatibility and provides a workaround for a newly discovered bug in FreeBayes.

The reference contigs now have a single N placed on both 5' and 3' ends of the contig. The current version of FreeBayes was ignoring some reads that mapped exactly to the 5' end of a reference contig, and adding the N padding fixes this issue.

Also, dDocent now outputs the reference in the correct orientation. Previous versions output reverse complemented reference contigs.  This had no impact on SNP calling; it's just important to note for comparisons to previous runs.

**It is recommended that all previous analyses are rerun with the current version.  The new version will most likely result in more SNP calls with higher coverage.**
