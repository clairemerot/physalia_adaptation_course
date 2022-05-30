# SNPs and SVs calling from whole genome resequencing data
Today we have explored the importance of structural variants in ecology and evolution, and you have tested different approaches to investigate haploblocks identified through reduced representation approaches such as RADseq. However, explicitely testing for the presence of SVs associated with haploblocks is not really possible with RADseq data (but see Yann's paper on detecting CNVs from RADseq data; Dorant et al. 2020, Molecular Ecology), so we'll now use whole genome resequencing data (Illumina short reads) from a few capelin samples. From these data, we can call both SNPs (sequence variation) and SVs (structural variation) and test whether the two provide concordant patterns of variation and differentiation, or not.

## Sequence variation
As both sequence and structural variant calling take quite some time, I did this step for you. Whether you want to try it running it when you have some spare time or use it in the future for the analysis of your own data, I will provide the command below. 
First, let's set up new folders for these analyses
```
cd
mkdir wgr #for whole genome resequencing
mkdir wgr/svs_delly
mkdir wgr/snps_bcftools
```

### 1. SNPs calling
Copy bam files (symbolic links) storing aligned reads to reference genome
```
cd wgr
ln -s /home/ubuntu/Share/WGS_bam/.bam* .
cp /home/ubuntu/Share/WGS_bam/bams_list.txt .
```
Call variants using bcftools
```
cd snps
bcftools mpileup -Ou -f ~/Share/resources/genome_mallotus_dummy.fasta -b bams_list.txt -q 5 -I -a AD,DP,SP,ADF,ADR -d 200 | bcftools call - -mv -Ov > snps_bcftools/capelin_wgs_unfiltered.vcf
```
Variant filtering is done in two steps. First we use 'bcftools filter' to filter SNPs based on mapping quality. Then, we apply a series of filters using VCFtools 
```
bcftools filter -e 'MQ < 30' snps_bcftools/capelin_wgs_unfiltered.vcf -Ov > snps_bcftools/capelin_wgs_filtered.tmp.vcf
### Count number of unfiltered SNPs surviving filters
grep -v ^\#\# snps_bcftools/capelin_wgs_filtered.tmp.vcf | wc -l  

echo "
>>> Filtering through VCFtools now!!
"
vcftools --vcf snps_bcftools/capelin_wgs_filtered.tmp.vcf \
    --minQ 30 \
    --minGQ 20 \
    --minDP 3 \
    --mac 2 \
    --max-alleles 2 \
    --max-missing 0.7 \
    --maf 0.05 \
    --recode \
    --stdout > snps_bcftools/capelin_wgs_filtered.vcf

### Count number of SNPs surviving filters
grep -v ^\#\# snps_bcftools/capelin_wgs_filtered.vcf | wc -l 
```

What is the proportion of SNPs surviving filters? How many are left? According to your dataset and needs, you may want to change or tweak these filters.

### 2. Analysis of sequence variation
As we have highlighted several times, a VCF is a VCF, whether the variants are called from RADseq or WGR data, whether you call hundreds or millions of variants, so you can repeat all the analyses we presented to you during the course using this new WGR dataset. For this tutorial we will focus on the visualization of genomic variation with PCA (see below).

## Structural variation
Calling structural variants (i.e. inversions, deletions, insertions, duplications, translocations) requires analyzing different types of information from the sequence data, such as overlap, orientation and read splitting. Therefore, we need software specifically designed to extract this kind of information. One commonly used of these programs is Delly2, which takes into account all of these 3 types of information. 


