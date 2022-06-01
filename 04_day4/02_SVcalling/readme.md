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
bcftools mpileup -Ou -f ~/Share/ressources/genome_mallotus_dummy.fasta -b bams_list_rg.txt -r Chr1 -q 5 -I -a AD,DP,SP,ADF,ADR -d 200 | bcftools call - -mv -Ov > snps_bcftools/capelin_wgs_unfiltered.vcf

```
Variant filtering is done in two steps. First we use 'bcftools filter' to filter SNPs based on mapping quality. Then, we apply a series of filters using VCFtools. Note that I'm calling variants from only chromosome 1 to save time. 
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
Also, because there was a problem with the internal individual flags in the read files, which cause the addition of an individual, I used this line to fix it
```
grep "##" capelin_wgs_filtered.vcf > capelin_wgs_filtered_fixed.vcf
grep -v "##" capelin_wgs_filtered.vcf | awk '{$10=""; print $0}' >> capelin_wgs_filtered_fixed.vcf
```

Because the variant calling takes quite some time, you can copy the original unfiltered and clean files to your working directory to explore the different files, theyir differences and play with filtering if you want
```
cp ~/Share/WGS_bam/snps_bcftools/capelin_wgs_*.vcf ~/wgr/snps_bcftools/.
```
What is the proportion of SNPs surviving filters? How many are left? According to your dataset and needs, you may want to change or tweak these filters.

Finally, we prepare the data in the same way as the SVs below and recode the genotypes in 012 format with 
```
vcftools --vcf capelin_wgs_filtered_fixed.vcf --012 --out capelin_wgs_filtered_fixed
```
and download the resulting files on your computer.

### 2. Analysis of sequence variation
As we have highlighted several times, a VCF is a VCF, whether the variants are called from RADseq or WGR data, whether you call hundreds or millions of variants, so you can repeat all the analyses we presented to you during the course using this new WGR dataset. For this tutorial we will focus on the visualization of genomic variation with PCA (see below).

## Structural variation
1. SVs calling
Calling structural variants (i.e. inversions, deletions, insertions, duplications, translocations) requires analyzing different types of information from the sequence data, such as read overlap, orientation and splitting. Therefore, we need software specifically designed to extract this kind of information. One commonly used of these programs is Delly2, which takes into account all of these 3 types of information. 
This below is the general approach to call SVs from Delly:
```
###Germline SV calling
#SV calling is done by sample for high-coverage genomes or in small batches for low-coverage genomes
delly call -g hg19.fa -o s1.bcf -x hg19.excl sample1.bam

#Merge SV sites into a unified site list
delly merge -o sites.bcf s1.bcf s2.bcf ... sN.bcf

#Genotype this merged SV site list across all samples. This can be run in parallel for each sample.
delly call -g hg19.fa -v sites.bcf -o s1.geno.bcf -x hg19.excl s1.bam

delly call -g hg19.fa -v sites.bcf -o sN.geno.bcf -x hg19.excl sN.bam

#Merge all genotyped samples to get a single VCF/BCF using bcftools merge
bcftools merge -m id -O b -o merged.bcf s1.geno.bcf s2.geno.bcf ... sN.geno.bcf

#Apply the germline SV filter which requires at least 20 unrelated samples
delly filter -f germline -o germline.bcf merged.bcf

```
However, in our case, we have low-coverage whole genome resequencing data from only 12 individuals, so we can call SVs directly from all the samples combine and forego the filtering step. Keep in mind that this is just a toy dataset, and the quality of these SV calls may not be very high.
```
cd
mkdir delly
cd delly
mkdir bams
ln -s /home/ubuntu/Share/WGS_bam/* .

### Run delly to call SVs on all samples combined
delly call -g ~/Share/resources/genome_mallotus_dummy.fasta -o capelin_sv.bcf BELB9_1.trimmed.sorted.bam BELD3_1.trimmed.sorted.bam BLA13_1.trimmed.sorted.bam BLA15_1.trimmed.sorted.bam BLA16_1.trimmed.sorted.bam BLA17_1.trimmed.sorted.bam BLA22_1.trimmed.sorted.bam BLA24_1.trimmed.sorted.bam BSO17_1.trimmed.sorted.bam BSO23_1.trimmed.sorted.bam BSO28_1.trimmed.sorted.bam POR19_1.trimmed.sorted.bam 

###Convert file from .bcf to .vcf
bcftools convert -O v -o capelin_sv.vcf capelin_sv.bcf
```
You can run this code, but it will take ~45 minutes. Because we don't have that time, you can copy the VCF file containing all the SVs into ~/wgr/svs_delly
```
cp ~/Share/WGS_bam/svs_delly/capelin_sv.vcf ~/wgr/svs_delly/.
```
2. VCF filtering and splitting
Though we don't have enough samples to run `delly filter` properly, we can do some filtering using the same approach we use for SNP filtering but in this case we'll filter just by missing data
```
vcftools --vcf capelin_sv.vcf \
    --max-missing 0.7 \
    --recode \
    --stdout > capelin_sv_filtered.vcf
```
Then, we can split the vcf file by SV type (or extract one particular SV type of interest with this script
```
grep "#" capelin_sv_filtered.vcf > capelin_sv_ins.vcf && grep "SVTYPE=INS" capelin_sv.vcf >> capelin_sv_ins.vcf ###insertions
grep "#" capelin_sv_filtered.vcf > capelin_sv_del.vcf && grep "SVTYPE=DEL" capelin_sv.vcf >> capelin_sv_del.vcf ###deletions
grep "#" capelin_sv_filtered.vcf > capelin_sv_inv.vcf && grep "SVTYPE=INV" capelin_sv.vcf >> capelin_sv_inv.vcf ###inversions
grep "#" capelin_sv_filtered.vcf > capelin_sv_dup.vcf && grep "SVTYPE=DUP" capelin_sv.vcf >> capelin_sv_dup.vcf ###duplications
grep "#" capelin_sv_filtered.vcf > capelin_sv_bnd.vcf && grep "SVTYPE=BND" capelin_sv.vcf >> capelin_sv_bnd.vcf ###break points
```
and count the number of SV identified, overall or for each type with
```
grep -v "#" capelin_sv_ins.vcf | wc -l
```
What is the most abundant SV type?

Next, we'll convert the genotypes in 012 format for PCA of the whole dataset and each specific 
```
###Recode genotypes for PCA
module load vcftools
vcftools --vcf capelin_sv_filtered.vcf --012 --out capelin_sv

###and for each SV type
vcftools --vcf capelin_sv_ins.vcf --012 --out capelin_sv_ins
vcftools --vcf capelin_sv_del.vcf --012 --out capelin_sv_del
vcftools --vcf capelin_sv_inv.vcf --012 --out capelin_sv_inv
vcftools --vcf capelin_sv_dup.vcf --012 --out capelin_sv_dup
vcftools --vcf capelin_sv_bnd.vcf --012 --out capelin_sv_bnd
```
and download the resulting files on your computer.


### PCA to compare patterns from different types of genetic variation
On your computer in R, perform one PCA and plot results based on each type of structural variant. The one below is the code for the PCA based on all SVs. Once you've done this, modify the script to perform the analysis on each type of SV and on the SNPs.
Note that the bam files that I used for thw two types of analysis, sequence and structural variation, were slightly different because of the introduction of an individual flag during mapping, which was later fixed. So, when you run the PCA with the SNPs dataset, you have to replace `Sample` with `Id`.
```
#load the population map with population assignment for each individual
popmap_delly<-read.table("popmap_capelin_wgr_delly.txt", header= TRUE)

#load geno data for SVs
geno_capelin_sv <- read.table("capelin_sv.012")[,-1] #load genotype matrix
geno_capelin_sv.pos <- read.table("capelin_sv.012.pos") %>% #load SNPs info
  mutate(., locus=paste(V1,V2,sep='_')) #create a new column for SNP info name (CHR + position)
geno_capelin_sv.indv <- read.table("capelin_sv.012.indv") #load individuals info

#Set rownames and colnames to the geno matrix
dimnames(geno_capelin_sv) <- list(geno_capelin_sv.indv$V1, geno_capelin_sv.pos$locus)
#check the geno matrix
geno_capelin_sv[1:12,1:10]

#impute missing data
geno_capelin_sv[geno_capelin_sv == -1] <- NA
geno_capelin_sv.imp <- apply(geno_capelin_sv,2,function(x){
  replace(x, is.na(x), as.numeric(names(which.max(table(x))))) })

geno_capelin_sv.imp[1:12,1:9]

##run and visualized PCA
pca.geno_capelin_sv <- prcomp(data.imp)
screeplot(pca.geno_capelin_sv)

#get stats info from the pca
sum.pca <- summary(pca.geno_capelin_sv)
#print stats info
sum.pca$importance[,1:5]

#prepare dataset to plot PCAs
pca.geno_capelin_sv.sub <- pca.geno_capelin_sv$x[,1:4] %>% #retain the first four PCs
  as.data.frame(.) %>% #transform to dataframe object
  tibble::rownames_to_column(., var='Sample') %>% #set rownames to a new column for samples ids
  dplyr::left_join(., popmap_delly, by='Sample')
pca.geno_capelin_sv.sub2<-cbind(pca.geno_capelin_sv.sub[1:5], popmap_delly[2:3])

#Here we use the left_join function
#from dplyr to wrap the population vector
#of our samples.
ggplot(pca.geno_capelin_sv.sub2) + aes(x=PC1, y=PC2, col=Pop) +
  ggtitle("PCA with all SVs from Delly") +
  geom_point() + 
  coord_fixed() +
  theme_bw() 
  ```

Do you see the same patterns from each type of variant? Let's discuss.
