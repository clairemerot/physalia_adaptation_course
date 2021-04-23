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
We will keep a list of the most outliers SNPs differentiating the haplotypes. Because BB is a small group, differentiation AA vs BB may highlight false outliers. We will thus rather extract the SNPs with the highest Fst between AA and AB.

```
AA_AB<-read.table("04_divergence/AA_AB.weir.fst", header=T)
head(AA_AB)
outlier_Chr4<-AA_AB[AA_AB$WEIR_AND_COCKERHAM_FST>=0.1,]
write.table(outlier_Chr4, "04_divergence/outlier_Chr4.txt", row.names=F, quote=F, sep="\t")
```
