# Physalia Adaptation Genomics Course

## Welcome to the Physalia Adaptation Genomics 2021 GitHub page
This repository includes scripts and data associated with the practical sessions of the Physalia Course about Adaptation Genomics
Prepared by Anna Tigano, Yann Dorant and Claire Mérot for the course in May 2021

Except for day 1 in which you won't have the bam files, all the tutorial can be completed with the compressed files stored in the github page for each day (day 2 to day 5). Each daily tutorial can be run independantly of the other days.

Here is the tentative schedule for the week. Within each block we'll keep some flexibility to have plenty of time for questions and discussion.
![schedule](schedule 2021.png)


## Logging on the AWS server from your computer

[AWS for Mac OS X and Linux users](https://github.com/clairemerot/physalia_adaptation_course/blob/2021/AWS_mac.md)


[AWS for windows users](https://github.com/clairemerot/physalia_adaptation_course/blob/2021/AWS_windows.md)


## Day1: Handling NGS data: From raw reads to SNPs matrix

Capelin: Data taken from Cayuela et al,2020. Molecular Ecology https://doi.org/10.1111/mec.15499
 
Genome assembly: For this course, we made a dummy assembly of about 90 MB (instead of about 500 MB) and 5 chromosomes (instead of 24).

Raw reads: Obtained by GBS methods (= Rad-seq), sequenced with IonTorrent.

1-1: Getting familiar with Unix environment

1-2: From raw sequences to mapped reads

1-3: Calling variants with Stacks
 
## Day2: Population structure and confounding factors

2-1: Fst statistics with vcftools (optionnal: with Stacks, optional: Pairwise-Fst and isolation-by-Distance )

2-2: Principal component analysis (PCA)

2-3: Discriminant Analysis of Principal Components (DAPC)

## Day3: Outlier detection - Environmental associations

Data: We focus on 12 population from Canada for which there is almost no geographic structure but environmental variability

3-1: Genetic structure and LD-pruned data

3-2: Outlier of differentiation with two methods (Outflank & Baypass)

3-3: Genotype-Environnement Associations with two methods (Baypass & redundancy analysis)

## Day4: Accounting for Structural Variants

We focus on 12 population from Canada. We recommend that you pick one of the two tutorials (haploblocks by local PCA or CNVs on RAD-seq data)

4-1: Investigating haplotypes blocks ( ~inversions?)

This tutorial include working on local PCA, but also calculation of LD, Fst and observed fraction of heterozygotes which may be useful in other contexts

4-2: Filtering duplicated loci in RAD-seq data ( ~ Copy number variants)

This tutorial show how to filter RAD loci to exclude duplicated ones (keep a reliable dataset for SNP analysis), and then how to analyse the duplicated loci for environmental associations.

## Day5: Functional approaches

5-1: SNPeff annotation of SNPs for coding & regulatory regions

5-2: Intersection between SNPs and genes with bedtools

5-3: Gene ontology enrichment

5-4: (Optional) Intersection between CNVs and repeats/TE
