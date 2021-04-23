## Step 4 Studying heterozygosity

We will now try to figure out whether heterozygosity is indeed higher in our Ab groups. We will use the --hardy options for vcftools which tests hardy-weinberg equilibrium for each SNP and report the observed and expected fraction of heterozygotes at that position

#### On the server: get H-W stats and number of heterozygotes
We will  use vcftools to create a vcf for each group 

```
vcftools --vcf 02_data/canada.vcf --keep 02_data/AA.list --recode --out 02_data/AA
vcftools --vcf 02_data/canada.vcf --keep 02_data/AB.list --recode --out 02_data/AB
vcftools --vcf 02_data/canada.vcf --keep 02_data/BB.list --recode --out 02_data/BB
ls 02_data
```

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
AB.hwe$het_fraction<-AB.hwe$Het_Obs/(AB.hwe$Homo1_obs+AB.hwe$Het_Obs+AB.hwe$Homo2_Obs)

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
