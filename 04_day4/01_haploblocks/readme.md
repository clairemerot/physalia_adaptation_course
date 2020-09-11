# Exploring the haploblock on Chr 4

## Prepare data for this tutorial

Thanks to our local PCA exploration on day 2, we know that there are non-recombining haploblocks which may be an inversion on chromosome 4.
We located the breakpoints approximately from 4.8MB to 16.6MB

Please copy all the folder 04_day4 from Share into you repository.

Please open the following path and we will run all commands from here
``` cd 04_day4/01_haploblocks```

You should have in the 02_data folder the vcf for canadian populations with all SNPs (1 snp randomly selected by locus) and all individuals.
You also have a file with information about individuals info_samples_canada.txt

``` ls 02_data```

## Back to the PCA: genotype the individuals

First, we want to classify our individuals according to the cluster that we observed on the PCA.
To do so, we want to extract the relevant regions, perform a PCA and clustering approach

We know how to exclude chromosome from a vcf. We will now select just a determined region from the vcf

#### on the server : make a vcf with chr 4 4.8-16.6MB
We will also export it as 012 as we did on day2 to do the PCA

```
gunzip 02_data/canada.vcf.gz
vcftools --vcf 02_data/canada.vcf --chr Chr4 --from-bp 4800000 --to-bp 16600000 --recode --out 02_data/canada.chr4inv
vcftools --vcf 02_data/canada.chr4inv.recode.vcf --012 --out 02_data/canada.chr4inv
```

Now you can copy the whole folder 01_haploblocks on your local computer in which we will run the PCA in R

#### on your computer in R studio

Like on day 2, we will perform a PCA (including earlier steps of filling missing data etc). Note that we could also have used the vcfR library as we did on day 3

```
geno <- read.table("02_data/canada.chr4inv.012")[,-1] #load geno
indv <- read.table("02_data/canada.chr4inv.012.indv") #load individuals info
rownames(geno)<-indv[,1]
geno[1:6,1:6] #check the geno matrix
geno.pca<-prcomp(geno) #run the pca
plot(geno.pca$x[,1],geno.pca$x[,2] # plot the pca
```
We tend to see three groups although the smallest one is not well define, possibly because this are less individuals belonging to this group.
To cluster those individuals into three groups, we will use kmeans methods

```
geno_kmean<-kmeans (geno.pca$x[,1], c(min(geno.pca$x[,1]), (min(geno.pca$x[,1])+max(geno.pca$x[,1]))/2, max(geno.pca$x[,1]) ))
geno_kmean
plot(geno.pca$x[,1],geno.pca$x[,2], col= geno_kmean$cluster, pch=20)
```
![pca](06_images/pca_cluster.png)

You can note the ratio of the between sum of squares over the total sum of square which is indicative of the fit of the clusters
Here it is around 92%, this is not amazing but still meaningful, clustering in 3 is better than a whole pool...

Let's keep those clusters as it is for subsequent analysis. In real life, you may want to refine further, for instance using windows or SNPs most strongly divergent, or excluding dubious intermediate individuals
To split the vcf between our three groups we will need a list of individuals id in each category

```
#save information
info_cluster_inv<-cbind(indv,geno_kmean$cluster)
colnames(info_cluster_inv)<-c("id_inv","cluster_inv")
head(info_cluster_inv)

#prepare the list
AA_ind<-info_cluster_inv[info_cluster_inv$cluster_inv=="1",1]
AB_ind<-info_cluster_inv[info_cluster_inv$cluster_inv=="2",1]
BB_ind<-info_cluster_inv[info_cluster_inv$cluster_inv=="3",1]
  
#export files
write.table(info_cluster_inv, "02_data/info_cluster_inv.txt", quote=F, row.names=F)
write.table(AA_ind, "02_data/AA.list", quote=F, row.names=F, col.names=F)
write.table(AB_ind, "02_data/AB.list", quote=F, row.names=F, col.names=F)
write.table(BB_ind, "02_data/BB.list", quote=F, row.names=F, col.names=F)
```

#### On the server
Please copy back those 4 files into the 01_haploblocks/02_data folder on the server for all subsequent analysis.
We will now use vcftools to create a vcf for each group (and simplify names)

```
vcftools --vcf 02_data/canada.vcf --keep 02_data/AA.list --recode --out 02_data/AA
vcftools --vcf 02_data/canada.vcf --keep 02_data/AB.list --recode --out 02_data/AB
vcftools --vcf 02_data/canada.vcf --keep 02_data/BB.list --recode --out 02_data/BB
ls 02_data
```

## Studying linkage disequilibrium

#### On the server: calculate Ld with Plink
To calculate LD we will use plink, but not in pruning mode. We want all pairwise LD on each chromosome.
for Ld, we will remove SNPs at low frequency as they will be uninformative and increase the size of the matrix (>5% of frequency - we could have filter up to 5% or 10% with whole genome data). 
We will focus on chromosome 4 but feel free to try other chromosome.

Plink requires three inputs (.bed, .bim, .fam). The argument --r2 calculate the Ld as R2 (you could also have chosen D), inter-chr makes a long matrix (square would have amke a suare one)
we need to put --allow-extra-chromosome since we are not on humans! and --ld-window-r2 0 asks all output to be printed. To reduce the file you can choose here a minimum threshold for R2

```
#extract a reduced vcf 
vcftools --vcf 02_data/canada.vcf --chr Chr4 --maf 0.05 --recode --out 03_ld/maf0.05_chr4

#format for plink
plink --vcf 03_ld/maf0.05_chr4.recode.vcf --make-bed --out 03_ld/maf0.05_chr4

#calculate ld
plink --bed 03_ld/maf0.05_chr4.bed \
--bim 03_ld/maf0.05_chr4.bim \
--fam 03_ld/maf0.05_chr4.fam \
--r2 inter-chr --allow-extra-chr --ld-window-r2 0 \
--out 03_ld/maf0.05_chr4

head 03_ld/maf0.05_chr4.ld
```

Let's do LD within a group of homokaryotes. We choose the group that has many individuals. For me that was AA.

We want to consider the same SNPs so we used the recode vcf >5% and keep the AA individuals using the command keep.

```
#extract a reduced vcf with the same snps but only BB individuals
vcftools --vcf 03_ld/maf0.05_chr4.recode.vcf --keep 02_data/AA.list --recode --out 03_ld/AA_maf0.05_chr4

#format for plink
plink --vcf 03_ld/AA_maf0.05_chr4.recode.vcf --make-bed --out 03_ld/AA_maf0.05_chr4

#calculate ld
plink --bed 03_ld/AA_maf0.05_chr4.bed \
--bim 03_ld/AA_maf0.05_chr4.bim \
--fam 03_ld/AA_maf0.05_chr4.fam \
--r2 inter-chr --allow-extra-chr --ld-window-r2 0 \
--out 03_ld/AA_maf0.05_chr4

head 03_ld/AA_maf0.05_chr4.ld
```


####On your computer : plotting LD

Please take the two files .ld on your local computer (in the 03_ld folder)
We will use R to visualise our results

```
library(ggplot2)
#load data
chr4.ld<-read.table("03_ld/maf0.05_chr4.ld", header=T)
head(chr4.ld)

#plot (very simple solution with ggplot. maybe you can find something nicer..?
ggplot(chr4.ld,aes(x=BP_A,y=BP_B, col=R2)) + theme_classic() + geom_point(size=1, shape=15) + 
  scale_colour_gradientn(colours=c("lightgrey","deepskyblue3","blue","blue3","navyblue","black"), limits=c(0,1),  name="R2")
 
#load data homozygotes
AA_chr4.ld<-read.table("03_ld/AA_maf0.05_chr4.ld", header=T)
head(AA_chr4.ld)

ggplot(AA_chr4.ld,aes(x=BP_A,y=BP_B, col=R2)) + theme_classic() + geom_point(size=1, shape=15) + 
  scale_colour_gradientn(colours=c("lightgrey","deepskyblue3","blue","blue3","navyblue","black"), limits=c(0,1),  name="R2")

#plotting both heatmap on the same graph...
ggplot(chr4.ld,aes(x=BP_A,y=BP_B, col=R2)) + theme_classic() + 
  geom_point( shape=15) + 
  geom_point(data=AA_chr4.ld, aes(x=BP_B,y=BP_A,col=R2), size=1, shape=15)+
  scale_colour_gradientn(colours=c("lightgrey","deepskyblue3","blue","blue3","navyblue","black"), limits=c(0,1),  name="R2")

```
![ld](06_images/Ld_heatmap.png)

What do you think? do you observe the linkage possibly due to an inversion (or a non-recombining block?)? Is it also observed in the BB group?

## Studying divergence (Fst)
to calculate Fst between our groups, we will do exactly as your did on day 2 with vcftools.
Note that here, this is not ideal since it is better to have balanced sample size (and our group AA is pretty small).

#### On the server: get Fst and FSt by windows
```
vcftools --vcf 02_data/canada.vcf --weir-fst-pop 02_data/AA.list --weir-fst-pop 02_data/AB.list --out 04_divergence/AA_AB
vcftools --vcf 02_data/canada.vcf --weir-fst-pop 02_data/AA.list --weir-fst-pop 02_data/BB.list --out 04_divergence/AA_BB
vcftools --vcf 02_data/canada.vcf --weir-fst-pop 02_data/AB.list --weir-fst-pop 02_data/BB.list --out 04_divergence/AB_BB
```

We can also do it by windows:
```
WINDOW=100000
WINDOW_STEP=25000
vcftools --vcf 02_data/canada.vcf --weir-fst-pop 02_data/AA.list --weir-fst-pop 02_data/AB.list --fst-window-size $WINDOW --fst-window-step $WINDOW_STEP --out 04_divergence/AA_AB
vcftools --vcf 02_data/canada.vcf --weir-fst-pop 02_data/AA.list --weir-fst-pop 02_data/BB.list --fst-window-size $WINDOW --fst-window-step $WINDOW_STEP --out 04_divergence/AA_BB
vcftools --vcf 02_data/canada.vcf --weir-fst-pop 02_data/AB.list --weir-fst-pop 02_data/BB.list --fst-window-size $WINDOW --fst-window-step $WINDOW_STEP --out 04_divergence/AB_BB
```

#### On your computer: visualise

```
#load data
AA_BB<-read.table("04_divergence/AA_BB.weir.fst", header=T)
head(AA_BB)

#plotFst_AAvsBB
ggplot(AA_BB, aes(x=POS/1000000, y=WEIR_AND_COCKERHAM_FST, col=CHROM))+
  geom_point()+ geom_smooth()+
  theme_classic()+
  facet_grid(cols = vars(CHROM), scales = "free_x", space="free_x")+
  labs(  x = "position (in MB)")
  
#load data
AA_BB.win<-read.table("04_divergence/AA_BB.windowed.weir.fst", header=T)
head(AA_BB.win)
```
![fst](06_images/Fst_AAvsBB.png)
As you can note within our region of interest on Chr4, some SNP have a super high Fst (up to 1), suggesting fixed alleles and extremely high divergence
You may have noticed that some FSt values are negatives.. This is likely driven by very low frequency alleles and inbalanced sample sizes. We also observe NA in the calculation of FSt by site

Try to plot also the AA_AB and AB_BB contrasts. 

We can also look at the values summarized by Fst
```
#calculate the position of the mid window
AA_BB.win$BIN_MID<-(AA_BB.win$BIN_START+AA_BB.win$BIN_END)/2

#plot
ggplot(AA_BB.win, aes(x=BIN_MID/1000000, y=MEAN_FST, col=CHROM))+
  geom_point()+ geom_smooth()+
  theme_classic()+
  facet_grid(cols = vars(CHROM), scales = "free_x", space="free_x")+
  labs(  x = "position (in MB)")
```

## Studying heterozygosity

We will now try to figure out whether heterozygosity is indeed higher in our Ab groups. We will use the --hardy options for vcftools which tests hardy-weinberg equilibrium for each SNP and report the observed and expected fraction of heterozygotes at that position

#### On the server: get H-W stats and number of heterozygotes
We will only keep SNPs with maf above 5% since rare SNP won't be super informative in terms of heterozygotes (high stochasticity due to sampling a low number of heterozygotes)

```
vcftools --vcf 02_data/BB.recode.vcf --maf 0.05 --hardy --out 05_heterozygosity/BB
vcftools --vcf 02_data/AB.recode.vcf --maf 0.05 --hardy --out 05_heterozygosity/AB
vcftools --vcf 02_data/AA.recode.vcf --maf 0.05 --hardy --out 05_heterozygosity/AA
```

You can look at the files ``` head 05_heterozygosity/AA.hwe ``` . As you see, it is really annoying to have a symbol "/" at the middle of our column. We will use a simple awk command to get a better format for subsequent R analysis

```
cat 05_heterozygosity/AA.hwe | awk -F"/" '$1=$1' OFS="\t" > 05_heterozygosity/AA_formatted.hwe
head 05_heterozygosity/AA_formatted.hwe
cat 05_heterozygosity/AB.hwe | awk -F"/" '$1=$1' OFS="\t" > 05_heterozygosity/AB_formatted.hwe
cat 05_heterozygosity/BB.hwe | awk -F"/" '$1=$1' OFS="\t" > 05_heterozygosity/BB_formatted.hwe
```

#### On your local computer: visualize results
Please copy the outputs (_formatted.hwe) on your local computer in the folder 05_heterozygosity and let's go into Rstudio to extract the information
For today, we will only be interested in the first 3 columns with the nnumber of observed homozygotes and heterozygotes but you may find useful to look at H-W stats for other purpose
So after loading the data, we will calculate the % of heterozygotes and plot that along the genome
And plot the data with our usual ggplot command

```
#load data
AB.hwe<-read.table("05_heterozygosity/AB_formatted.hwe", header = T)
#rename columns
colnames(AB.hwe)<-c("CHR","POS","Homo1_obs", "Het_Obs", "Homo2_Obs", "Homo1_Exp", "Het_Exp","Homo2_Exp","Chisq_HWE","P_HWE","P_HET_DEFICIT", "P_HET_EXCESS")
head(AB.hwe)
#calculate the fraction of observed heterozygotes
AB.hwe$het_fraction<-BB.hwe$Het_Obs/(BB.hwe$Homo1_obs+BB.hwe$Het_Obs+BB.hwe$Homo2_Obs)

#plot
ggplot(AB.hwe, aes(x=POS/1000000, y=het_fraction, col=CHR))+
  geom_point(alpha=0.5)+ geom_smooth()+
  theme_classic()+
  facet_grid(cols = vars(CHR), scales = "free_x", space="free_x")+
  labs(  x = "position (in MB)")
```
You can do the same for AA and BB. Now let's try to visualise the three of them together

```
#keep genotype information before joining the three tables
BB.hwe$geno<-"BB"
AB.hwe$geno<-"AB"
AA.hwe$geno<-"AA"
all.hwe<-rbind(AA.hwe, AB.hwe, BB.hwe)
head(all.hwe)

#plot, colouring by genotype
ggplot(all.hwe, aes(x=POS/1000000, y=het_fraction, group=geno, col=geno))+
  geom_point(alpha=0.5)+geom_smooth()+
  theme_classic()+
  facet_grid(cols = vars(CHR), scales = "free_x", space="free_x")+
  labs(  x = "position (in MB)")
```
![Hobs_all](06_images/Hobs_all.png)

We could also do violin plots by chromosome and by groups. And if you want to do something more advanced, you can try a violinplot of the region inside the inversion vs. outside 
```
ggplot(all.hwe, aes(x=geno, y=het_fraction, col=geno))+
  geom_violin()+
  stat_summary(fun=median, geom="point", size=2) +
  facet_grid(cols = vars(CHR))+
  theme_classic()
```
![Hobs_violin](06_images/Hobs_violin.png)

As you observed on the Manhattan plots, there is a lot of heterogeneity between SNPs. Perhaps it might be worth looking at results by sliding-windows?
Our case is not ideal because SNPs are sparese (RAD-seq) but with shole-genome data you would have no choice bu doing windows.

We will use an easy function provided in the package windowscannr. It is not super fast but with our data, that should be ok.
If you have not installed it already you can do with 
```
library(devtools)
install_github('tavareshugo/windowscanr')
library(windowscannr)
```
As argument you can give the window size, step and whether it is doing the mean (or another summary statistics) within each window.
We need to give groups = chromosomes to avoid joining positions on different chromosomes!

```
WINDOW<-500000
WINDOW_STEP<-100000
BB.hwe_win <- winScan(x = BB.hwe,groups = "CHR", position = "POS",values = c("het_fraction"),win_size = WINDOW,win_step = WINDOW_STEP,funs = c("mean"))
head(BB.hwe_win)
```

You can perform the window summary on AA and AB, try plotting as above group by group, or all groups together.
It is slightly more readable. 

In all cases, we nevertheless note the expected higher observed heterozygosity in AB around the middle of Chr4. 
On Chr5, there is a region of high heterozygosity in all three groups, which may be driven by sex.
