# Exploring heterogeneity in the genome and putative rearrangement on Chr 4
IMPORTANT: Please copy all the folder 04_day4 from /Share into you repository on the AWS server. Please copy it in your own computer and try to follow the same architecture when moving your files.

``` 
cp -r ~/Share/physalia_adaptation/04_day4 .
cd 04_day4/01_haploblocks
```
We will run all commands from here

## Step 0 Sliding PCA 
As you saw on day 2, the PCA performed on the 240 samples from North America display a very unexpected pattern. The loadings indicate that some portion of the genome are overwhelmingly driving the structure.We suspect there may be sex-linked markers and/or chromosomal rearrangements
To get a better sense of what's going on, we will be doing PCA again, but along the genome by window of X SNPs.For this we will use  R package * lostruct* available here https://github.com/petrelharp/local_pca and presented in this publication https://www.genetics.org/content/211/1/289

This analysis is more powerful if we keep all SNPs including those in LD so we will come back to the vcf unfiltered

####  prepare the vcf (in the terminal) [do not run]
We will jump over the preparation of the file because we did not manage to make R communicate with bcftools on the AWS for a mysterious reason. It works very well on our university cluster though... So keep in mind those preparative steps if you want to re-do the analysis on your dataset and for now, please read through it and start applying the tutorial only at the end of the jump.

####[EDIT JUMP -> don't run this part, please read it for your information and run the commands after END OF JUMP]####
First you need to process the unfiltered vcf for the 12 NWA populations.We will tend to use here the vcf before selecting one random SNP per loci (i. e. a vcf with all SNPs even if they are linked) because we need to be as dense as possible and we are looking for LD. I put a version of the full vcf in Share/ressources.
We will use several programms to convert the vcf into a bcf, sort it and index it.
The package lostruct can also take a vcf but it is heavy to load the whole file. Instead they made function that is able to cut the vcf window by window to avoid overloading the memory (but this requires a sorted, indexed bcf)
```
#we need to sort the vcf
(grep ^"#" ../../Share/ressources/capelin_NWA.vcf; grep -v ^"#" capelin_NWA.vcf | sort -k1,1 -k2,2n) > 00_localPCA/capelin_NWA_sorted.vcf

#then compress it and index it
bgzip -c 00_localPCA/capelin_NWA_sorted.vcf > 00_localPCA/capelin_NWA_sorted.vcf.gz
tabix -fp vcf 00_localPCA/capelin_NWA_sorted.vcf.gz

#then convert it to a bcf and sort it
bcftools convert -O b 00_localPCA/capelin_NWA_sorted.vcf.gz > 00_localPCA/capelin_NWA_sorted.bcf
bcftools index 00_localPCA/capelin_NWA_sorted.bcf
```
To check whether this has worked and produced files, you can do a quick ```ls -lh``` which will show you the size of the files in human-readable format

Now we are good to work in R with the library. This requires some computational power and memory, so we suggest to make the initial steps in R command lines on the server and then copy the output files to visualise on your local Rstudio

####  run lostruct (in the terminal)[do not run]
To start R in command line, just type "R". Now you have a R console and we wil run the lostruct procedure
```
#open library
library(lostruct)
options(datatable.fread.input.cmd.message=FALSE)  #disable a useless message
snps <- vcf_windower("00_localPCA/capelin_NWA_sorted.bcf",size=100,type='snp', sites= vcf_positions("00_localPCA/capelin_NWA_sorted.bcf"))
```
This function makes windows out of the given data file of your chosen size. You can choose the size of the window with "size" and on which variable you want to split ('snp' or 'bp'). We suggest to use window of 100 snp since we are not very dense (RAD data) and we don't have a lot of snps.Typically with whole genome you may first run by windows of 1000 or 5000 snps for a first look, and then refine with smaller windows. The analysis can be run chromosome by chromosome (as in the paper) or on the entire genome. Here, we are going for the entire genome.
You can display for instance the 5th window and know its location by doing
```
snps(5)
region(snps) (5)
```
We can now run the PCA on all windows. Here we choose to consider k=npc=2 because they usually capture most variance for each local PCA
```
pcs <- eigen_windows(snps,k=2)
dim(pcs) #check dimension
head (pcs[,1:10]) #look at the first 10 columns
pcs_noNA<-pcs[-which(is.na(pcs[,1])),] #because of NA, some windows were not computed by pca. we will remove them
```
In the matrix pcs, each rows give the first k eigenvalues and k eigenvectors for each window. This gives you a matrix with 483 columns (3 columns of info, 240 columns with PC1 score for each individual, and 240 column with PC2 score for each individual). It has as many rows as windows (1016 with windows of 100 SNPs)

As you see we don't know the position of each window, we will get it with the function regions, remove the NA windows and exapnd the pca matrix to include the position information we retrieve before and export the file
```
#retrieve positions
window_pos<-region(snps)()
head(window_pos)

#keep windows without NA
window_pos_noNA<- window_pos[-which(is.na(pcs[,1])),]
#merge
pca_matrix_noNA<-cbind(window_pos_noNA, pcs_noNA)
head (pca_matrix_noNA[,1:10])

#save the file
write.table(pca_matrix_noNA, "00_localPCA/pca_matrix.txt", sep="\t", row.names=FALSE, quote=FALSE)
```
####[END OF JUMP]#### Retrieve the file I did for you:

####  run lostruct (in the terminal)[yes run it!]

Back to work, now that you have read and understood, let's comme back to our terminal. I have put for you the pca_matrix.txt in the folder 05_localPCA. please load it into R, and split it into the window information and the PCs.To start R in command line, just type "R". Now you have a R console and we wil run the end of lostruct procedure
```
library(lostruct)
pca_matrix_noNA<-read.table("00_localPCA/pca_matrix.txt", sep="\t", header=T, stringsAsFactors=FALSE)
head(pca_matrix_noNA)
window_pos_noNA<-pca_matrix_noNA[,1:3]
pcs_noNA<-as.matrix(pca_matrix_noNA[,4:dim(pca_matrix_noNA)[2]])
```

Now we can go back to the main pipeline:

The lostruct procedure proposes to compute pairwise distances between those windows and visualise it with a MDS (multidimensional scaling). Our goal is to identify groups of windows which display similar PCA pattern.This is done with the following functions (we uses 2 PC per window as above, and will look at the 1st 10 axes of the MDS)
```
pcdist <- pc_dist(pcs_noNA,npc=2)
mds_axe<-cmdscale(pcdist, k=10)
head(mds_axe)

#again the mds file is missing position information so:
mds_matrix<-cbind(window_pos_noNA, mds_axe)
write.table(mds_matrix, "00_localPCA/mds_matrix.txt", sep="\t", row.names=FALSE, quote=FALSE)
```
Since this is a little long , we will let it run and explore on your local computer with Rstudio some of the local pca. Let's download pca_matrix.txt locally and we will play in R studio to look at some of those local PCA

#### Visualising the local PCA outputs (on your computer in Rstudio)
Back on our local computer in R studio, we will look at all those local PCA.
Set your working directory as 04_day4/01_haploblocks, load useful libraries (ggplot2) and the matrix of pca

```
pca_matrix<-read.table("00_localPCA/pca_matrix.txt", header=TRUE)
pca_matrix[1:10,1:10]
n_windows<-dim(pca_matrix)[1] # the number of windows we have
```

We may want to simply plot the pca for some windows. This is not the most easy because remember the format is a bit tricky
Look at the format. We have 3 columns for position, total eigen values, eigvalue of PC1, of PC2 and then 240 values for PC1 scores of all our samples, and 240 values for PC2 scores of all samples
```
   chrom   start     end      total    lam_1     lam_2   PC_1_L_01    PC_1_L_02    PC_1_L_03    PC_1_L_04
1   Chr1    4598  102627 11.2676700 2.237855 1.2467539  0.04757381  0.003211281  0.048236420 -0.029799279
2   Chr1  102627  298945  5.4675252 1.664517 1.1009505  0.05177397 -0.066433281  0.056306263  0.057187867
3   Chr1  298949  557448 20.3195707 3.514855 2.5803417 -0.07997331 -0.014285140 -0.077618567  0.061921540
4   Chr1  557463  664701  7.9826780 1.707709 1.4688195  0.00317093 -0.027442327 -0.003586276 -0.025939893
```

So to get information for a window we can fo something like:
```
Nind<-240

pc1_i<-t(pca_matrix[i, 7:(Nind+6)]) #scores along PC1
pc2_i<-t(pca_matrix[i, (Nind+7):(2*Nind+6)]) #scores along PC2
var1<-round(pca_matrix[i, 5]/pca_matrix[i, 4],2)*100 # % of variance explained by PC1
var2<-round(pca_matrix[i, 6]/pca_matrix[i, 4],2)*100 # % of variance explained by PC2
midpos_i<-(pca_matrix[i, 2]+pca_matrix[i, 3])/2 #average position of the window
window_i<-paste(pca_matrix[i, 1], midpos_i , sep="_") #paste the name of CHR and the midposition

plot(pc1_i, pc2_i, pch=20, xlab=paste("PC1", var1 , "%"), ylab=paste("PC2", var2, "%"), main=window_i)
```

Not super clean but it works. Building on  that you can do anything to reformat your matrix of PC, eigen values, etc...

Now what do we want to know, we want to look which windows explain the pattern observed in the global pca. We can look at correlation between global PCs and PC1 of each local PCA. I propose to take the PCA performed on the 12 NWA populations, and look at the correlation between PC1 of the global PCA and PC1 of each local PCA. (then you can do the same with PC2 of the global PCA and PC1 of each local PCA...) To help, I put the geno-transformed matrix inside your folder so that you don't go back to vcftools but the .012 file has been done exactly as shown on day2

I suggest below a very basic loop to store the correlation by windows. You can probably do something more fancy :-)
Please note here that we could also have used the genotype correlation, for instance call 0/1 or 0/1/2 the cluster observed on the global pca and then look for each snps at the correlation between genotypes and gneotype for the cluster identified.

```
geno <- read.table("canada.012")[,-1] #load geno
geno[1:6,1:6] #check the geno matrix
global.pca<-prcomp(geno) #run the pca
plot(geno.pca$x[,1],geno.pca$x[,2] # plot the pca
PC_of_interest<-global_pca$x[,1] #if you want to look at correlation with PC1

#initialise the vector
corr_vector<- vector(length=n_windows)
#loop over windows to store correlation factor
for (i in 1 : n_windows)
{
  pc1_i<-t(pca_matrix[i, 7:(Nind+6)]) #scores along PC1
  corr_vector[i]<-abs(cor(PC_of_interest, pc1_i)[1,1])
}

```
Now let's merge correlation and position to do a Manhattan plot of correlation along the genome. We need to take a midposition for each window.
I suggest that we use ggplot to visualize and facet_grid is a useful way to make quick Manhattan plot with chromosome side by side
```
pca_correlation<-cbind(pca_matrix[,1:3], corr_vector)
pca_correlation$midpos<-(pca_correlation$start+pca_correlation$end)/2
head(pca_correlation)

ggplot(pca_correlation, aes(x=midpos, y=corr_vector, colour=chrom))+
  geom_point()+
  theme_classic()+
  facet_grid(cols = vars(chrom), scales = "free_x", space="free_x")
```

What do you see?  Which windows correlate with PC1? What do you think?
![local_pca_corr_PC1](https://github.com/clairemerot/physalia_adaptation_course/blob/master/02_day2/05_localPCA/images/local_pca_corr_PC1.png)

Now you can have a look at correlation between local PC1s and the global PC2...

####  Using the MDS
Here that was easy, because we knew there was something weird on PC1 and PC2. But please keep in mind that, even if on the global PCA, no region is driving a specific clustering, there may still be, on some chromosome, or some regions, similar clustering of individuals that reveal population structure, chromosomal rearragements, sex, non recombining haploblocks, etc. Exploring the MDS is a way to detect such heterogeneity in the genome.

Let's load the mds and plot the first axes
```
mds_matrix<-read.table("mds_matrix.txt", header=TRUE)
head(mds_matrix)
mds_matrix$midpos<-(mds_matrix$start+mds_matrix$end)/2

#for ggplot you need to rename the columns
colnames(mds_matrix)<-c("chrom","start","end","mds1","mds2","mds3","mds4","mds5","mds6","mds7","mds8","mds9","mds10")

ggplot(mds_matrix, aes(x=mds1, y=mds2, colour=chrom))+
  geom_point()+
  theme_classic()
ggplot(mds_matrix, aes(x=mds3, y=mds4, colour=chrom))+
  geom_point()+
  theme_classic()
```
![mds1_2](https://github.com/clairemerot/physalia_adaptation_course/blob/master/02_day2/05_localPCA/images/mds1_2.png)

As you see, MDS 1 and 2 are largely driven by Chr 4 and Chr5. Let's look at mds scores along the genome to pinpoint those regions
Building on what we did before, you can probably make a Manhattan plot with a midposition as x, and mds1 or mds2 as y.

This is the output for MDS2
![mds2_manhattan_plot](https://github.com/clairemerot/physalia_adaptation_course/blob/master/02_day2/05_localPCA/images/mds2_position.png)

To follow-up, You can try to find approximately the breakpoints of those areas that appear outliers.
We can locate them approximately from 4.8MB to 16.6MB on chromosome 4 and the full chromosome 5

#### A note about running R into a terminal or on a server
Instead of running R frontally, as we did at the beginning for lostruct, we could have written a whole script and run it with
```
Rscript my_fancy_script.R "option1" "option2"
```
The first lines my_fancy_script.R would be:
```
argv <- commandArgs(T)
option1 <- argv[1]
option2 <- argv[2]
```
## Step 1 Genotype the individuals for the haploblocks
Thanks to our local PCA exploration on day 2, we know that there are non-recombining haploblocks which may be an inversion on chromosome 4.
We located the breakpoints approximately from 4.8MB to 16.6MB

First, we want to classify our individuals according to the cluster that we observed on the PCA.
To do so, we want to extract the relevant regions, perform a PCA and clustering approach

We know how to exclude chromosome from a vcf. We will now select just a determined region from the vcf

#### on the server : make a vcf with chr 4 4.8-16.6MB
You should have in the 02_data folder the vcf for canadian populations with all SNPs (1 snp randomly selected by locus) and all individuals.
You also have a file with information about individuals info_samples_canada.txt
``` ls 02_data```

We will make a vcf with chr 4 4.8-16.6MB and export it as 012 as we did on day2 to do the PCA

```
gunzip 02_data/canada.vcf.gz
vcftools --vcf 02_data/canada.vcf --chr Chr4 --from-bp 4800000 --to-bp 16600000 --recode --out 02_data/canada.chr4inv
vcftools --vcf 02_data/canada.chr4inv.recode.vcf --012 --out 02_data/canada.chr4inv
```

Now you can copy the generated files into  on your local computer following the same architecture (04_day4/01_haploblocks/02_data) 

#### on your computer in R studio
 perform a PCA 
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

## Step 2 Study linkage disequilibrium

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

## Step 3 Studying divergence (Fst)
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
We will keep a list of the most outliers SNPs differentiating the haplotypes. Because BB is a small group, differentiation AA vs BB may highlight false outliers. We will thus rather extract the SNPs with the highest Fst between AA and AB.

```
AA_AB<-read.table("04_divergence/AA_AB.weir.fst", header=T)
head(AA_AB)
outlier_Chr4<-AA_AB[AA_AB$WEIR_AND_COCKERHAM_FST>=0.1,]
write.table(outlier_Chr4, "04_divergence/outlier_Chr4.txt", row.names=F, quote=F, sep="\t")
```

## Step 4 Studying heterozygosity

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

If you are curious and if you have time, you may want to try exploring the Chr5 in the same way. Note that you won't need the first step to find clusters since we already have the sex information.
