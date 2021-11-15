# Karyinon (working title)

This workflow began as a longread version of the Karyon pipeline (below), and has now evolved into a full assembly analysis suite, designed particularly for analysis of complex polyploid genomes.


This workflow includes full assembly appraisal with the following methods:
-Read mapping with minimap2 and Vulcan
-QUAST assembly appraisal statistics
-mitochondrial contig flagging/removal
-contaminant identification with Blobtools
-homeologous contig detection using Blast, Nucmer, and Mummer
-Statistics on homeolog differences with DnaDiff
-Homeologous and haplotig purging with redundans
-Variant calling with GATK
-Flagstat generation with Samtools
-Ploidy estimations with Karyon workflow
-Coverage analytics and plots with Mosdepth

---
In the future it will also include:

-Merqury

---

Karyon is a package that assesses the ploidy of a library, reads or contigs.

Github repo: 

[Gabaldonlab/karyon](https://github.com/Gabaldonlab/karyon)

Paper: 

[Karyon: a computational framework for the diagnosis of hybrids, aneuploids, and other non-standard architectures in genome assemblies](https://www.biorxiv.org/content/10.1101/2021.05.23.445324v1?rss=1)
