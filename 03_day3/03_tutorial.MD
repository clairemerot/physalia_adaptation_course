# Day3: Detecting putatively adaptive loci
Today, we will restrain our analysis on the Canadian 12 populations which belong to the same lineage and for which we have environmental data.

To start, please copy the folder Day 3 into your home directory. We have put inside the folder 02_data your two vcf from yesterday:
- the vcf of the 12 populations filtered at 1 SNP per locus
- the vcf of the 12 populations filtered at 1 SNP per locus without chr 5 (sex linked) and without the putative rearranged regions on chr4
- an info file including sex, as well as population, latitude, longitude and temperature for the resticted dataset of 240 samples
- an info file about the 12 populations

Those files were all generated (by you!) on the previous days but to make you start more easily, we decided to prepare them here.
You can have a look at the files with the command 'head 01-data/file' (or less 02_data/file)
vcf files will need to be unzipped wiht the 'gunzip' command.
To see the beginning of the vcf file it may be useful to use 'less -S 02_data/file', which cut the right part of the file to fit your window, or 'head -n 25 02_data/file' to look at the first 25 lines

We will work from the 03_day3 folder and try to open our data within their folder (for instance 'head 02_data/file') and ouput in the folder of the analysis for instance (> 03_fst/output)

## 3A Investigate F-statistics in relation with geography
### Pairwise differentiation between populations

**[On the Amazon server]**
We will calculate here FST between all pairs of populations. FST accross the genome are expected to be largely driven by neutral markers while peaks of FSt may be related to selection for local adaptation... So we will first try to get a sense of global pattern and then look for outliers

To get pairwise FSt, we will use the R package STAMPP.
to save time today, we will use a toolbox developped by Yann Dorant. You may be interested in looking at the scripts to understand how this is done.
This toolbox embed various useful script in order to fastly convert vcf format to common pop genomics formats (genepop, StAMPP, baypass, bayenv...). If you are interested to learn more about this toolbox, you will find the full description at https://gitlab.com/YDorant/Toolbox

To download the toolbox in your current working folder on the server (03_day3), use the following command line:
```
git clone https://gitlab.com/YDorant/Toolbox
```

Ok, now we are ready to convert our VCF files to the StAMPP file format. The toolbox have an easy way to do that with a bash script. This bash script require four args:
* -v VCF file
* -p population map
* -f output file format
* -o output prefix name

```
bash Toolbox/00-VCF_Reshaper.sh -v 02_data/population.can.random.snp.vcf -p 02_data/popmap_canada.txt -f StAMPP -o 03_fst/population_can_random.snp
```
Then move inside the toolbox folder ``cd Toolbox/``
Once is done, you will ses that a ``.StAMPP`` file inthe 03_fst folder.
Then, we can run the Rscript `StAMPP-fst.R` to perform pairwise FST for each dataset. This script require three args.
* Input StAMPP file
* Output prefix
* Number of CPU allowed (default CPU=1)

```
Rscript Toolbox/StAMPP-fst.R 03_fst/population_can_random.snp.StAMPP 03_fst/population_can_random.snp 3
```
Allowing 3 CPU, each FST calculation should take around 5-6 minutes

Once FST calculations are done, you will see that four FST output files have been generated per dataset.
* prefix_fst_bootstrapes.txt
* prefix_fst_matrix.txt
* prefix_fst_pvalue.txt
* prefix_fst_reshape.txt

You can explore each of one with the cmd line ``less -S file.txt``. However, today, we only use the FST matrix.

So, now you can export the pairwise FST matrix (file suffix ``_fst_matrix.txt``) on your local computer to make a few visualisations in R
Please also export the info files about populations too
* 02_data/info_pop_geo_eco.txt
* 02_data/info_samples_canada.txt

**[On your local computer]**
We are now in Rstudio on your computer.
We keep it simple and do a simple numeric matrix but you can imagine more fancy ways, with heatmaps or so.

First, load the required libraries
```
# Libraries
  library(dplyr)
  library(magrittr)
  library(tibble)
  library(gplots)
  library(RColorBrewer)
  library(corrplot)
# --------------
```
Second, I give you a useful fonction to deal with the FST matrix. In fact, this is a triangular matrix and we need to fill the upper diag to use the heatmap function.
So, to do this, I give you the following function. Add it in your Rscript.
```
makeSymm <- function(m, position) {
    # Add symetrical triangle matrix (upper or lower)
  if (position == 'upper'){
    m[upper.tri(m)] <- t(m)[upper.tri(m)]
    return(m)
  }
  if (position == 'lower'){
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
    return(m)
  }
}
```
Now we can load the data, arrange it and then plot.

```
#-------------- fst matrix for all SNPs --------------------
fst.mat <- read.table("03_fst/population_can_random.snp_fst_matrix.txt")
#use the function to fill the full matrix
fst.all.mat<- fst.mat %>%
              as.matrix(.) %>%
              makeSymm(., 'upper')
fst.all.mat[is.na(fst.all.mat)] <-  0 #replace NAs by 0 (NAs unaccepted for the heatmap function)
fst.all.mat[1:10,1:10] #check the fst_matrix
              
#visualise values
corrplot(fst.all.mat, is.corr = FALSE, method="number", addgrid.col = FALSE, diag=F, type="lower", number.digits = 3, number.cex=0.7)

#Visualize pairwise FST through a heatmap
gplots::heatmap.2(fst.all.mat, trace = 'none',
                  col= colorRampPalette(brewer.pal(9, "Reds"))(15),
                  key.xlab='FST')
```
Note that values goes up to FST = 0.02 ! but most of them are very low
![img_fst_all](07_img_readme/Fst_heatmap_all_snps.png)

What do you notice? Is it heterogeneous? Do some population look more differentiated than other?
Why do you think A and J are so different? 

### Isolation by distance
To explore whether this is linked to distance between population, we will do a IBD test (Isolation by distance)
```
library(SoDA)
library(reshape2)
library(dplyr)
library(magrittr)
library(tibble)
library(ggplot2)

#import information about populations
info_pop <- read.table("02_data/info_pop_geo_eco.txt", header=T)
head(info_pop)
#calculate geogrpahic (euclidian) distances between all pairs of populations
distance <- dist(SoDA::geoXY(info_pop$latitude, info_pop$longitude)) %>%
  as.matrix(.)
dimnames(distance) <- list(info_pop$pop,info_pop$pop)
distance

#prepare datasets
#linearize the distance matrix
dist.melt <- reshape2::melt(distance) %>%
  set_colnames(., c('pop1', 'pop2','distance'))
head(dist.melt)

#linearize the fst matrix
fst.melt <- reshape2::melt(fst.all.mat) %>%
  set_colnames(., c('pop1', 'pop2','FST'))

#join the distance and fst
IBD.df <- left_join(dist.melt, fst.melt, by=c('pop1','pop2')) %>%
  filter(., distance > 0)
head(IBD.df)

#test association with FST
cor.test(log(IBD.df$distance), IBD.df$FST/(1-IBD.df$FST))

#plot IBD
ggplot(IBD.df) + aes(x=log(distance), y=FST/(1-FST)) +
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +
  theme_bw()
```
![img_IBD](07_img_readme/IBD_plot_all_snps.png)

So Isolation by distance is not significant and it does not seem that geography can explain the genetic distances very well with this full dataset.
If we come back to our heatmap, we can notice the cluster of populations C, F, I. How do you interpret it?

Let's look at the sex ratio in our data.
```
info_ind<-read.table("02_data/info_samples_canada.txt", header=T)
head(info_ind)
table (info_ind$pop,info_ind$sex)
```
```
    F  M
  A  0 20
  B 10 10
  C  0 20
  D 15  5
  E 10 10
  F  5 15
  G 10 10
  H 10 10
  I  0 20
  J 18  2
  K 10 10
  L 10 10
 ```
It seems that the field sampling has not been very good at balancing sex-ratio between populations... We should be worried about sex-linked markers driving the signal!It may be one of the reason why A and J are the most diferrentiated
(note: we intentionnaly subsetted the dataset to create this bias which was not in the original publication ;-) but it may happen easily for some species or low sample size)

### Getting rid of LD spurious effects?
Yesterday, we saw that a region of chromosome 4 and sex-linked markers on chr 5 were overwhelming the structure. Will that influence our pairwise estimates? Possibly
Let's re-run the steps above on the vcf in which we removed the chr5 and chr4.

Let's look at the pairwise FST matrix and the IBD stats
![img_fst_all](04_img_readme/Fst_heatmap_no_chr4-5.png)

What do you see now?
pay attention to the absolute value of FST, to the clustering and IBD (or the absence of)?

### Making a LD-pruned vcf
Another way to get rid of SNPs in LD due to the inversion or sex (?) would be to use a LD-pruned set of SNPs. For this we can use PLINK which will compute LD between SNPs by windows along the genome and keep one SNP out of several SNPs in linkage. Here we are interested in long distance LD (since we have already removed short-distance LD by keepin only one SNP per RAD locus) so I have set the window quite large.

Plink requires a bim file when taking a vcf input this can be created with the first line.
Then we give our parameters. We chose to be very stringent removing SNP with a VIF>2 in windows of 100 SNPs or 100KB
We could also have given a R² threshold with --indep-pairwise $WINDOW $SNP $R2

```
plink --vcf 02_data/population.can.random.snp.vcf.gz --make-bed --out 02_data/population.can.random.snp
 
WINDOW=100
SNP=100
VIF=2
plink --bed 02_data/population.can.random.snp.bed \
--bim 02_data/population.can.random.snp.bim \
--fam 02_data/population.can.random.snp.fam \
--indep $WINDOW['kb'] $SNP $VIF --allow-extra-chr \
--out 02_data/population.can.random.snp
```
You can have a look at the files 02_data/population.can.random.snp.prune.in and 02_data/population.can.random.snp.prune.outThey include a lit of SNP id as displayed in the vcf ("53:2:+") for instance. This list will be useful forour subsequent analyses in OutFLANK. For other purpose we may also want to have a vcf with only pruned marker. We are using again vcftools with the --exclude option (and recode)
```
vcftools --vcf 02_data/population.can.random.snp.vcf --exclude 02_data/population.can.random.snp.prune.out --recode --out 02_data/population.can.random.snp.pruned
```
If you have finished early, you can play at doing the pairwise FST/IBD or the population structure analysis of yesterday on the pruned dataset to check whether this removed the structure linked to sex or Chr4 rearrangement...


## 3B Investigate outliers of differentiation
### With OutFLANK
 OutFLANK is an R package that implements the method developed by Whitlock and Lotterhos (https://www.journals.uchicago.edu/doi/10.1086/682949) to use likelihood on a trimmed distribution of FST values to infer the distribution of FST for neutral markers. This distribution is then used to assign q-values to each locus to detect outliers that may be due to spatially heterogeneous selection.
https://github.com/whitlock/OutFLANK

Whitlock, M. C., and K. J. Lotterhos. 2015. Reliable detection of loci responsible for local adaptation: Inference of a neutral model through trimming the distribution of FST. The American Naturalist. 186:S24–S36.

It has a super good vignette here: https://htmlpreview.github.io/?https://github.com/whitlock/OutFLANK/blob/master/inst/doc/OutFLANKAnalysis.html
We will more or less follow it today.

#### Step 0 prepare the data
Let's first convert our vcf to the outflank format:
```
#we use the library vcfR to convert the vcf into the OutFLANK format
library(OutFLANK)
library(vcfR)
lirabry(ggplot2)
obj.vcfR <- read.vcfR("02_data/population.can.random.snp.vcf")

#extract useful informations about snp id and position
position <- getPOS(obj.vcfR) # Positions in bp
chromosome <- getCHROM(obj.vcfR) # Chromosome information
id_snp <- getID(obj.vcfR) # ID of the SNP

#gather this info in a dataframe
chr_pos<-as.data.frame(cbind(id_snp, chromosome, position)) # save info about id, chr, position
#R is sometimes not good at categorizing columns and here i had a problem that position was a factor... 
#this is an easy way to transform into a numeric
chr_pos$position<-as.numeric(as.character(chr_pos$position)) 

#we expect that it will be useful for subsequent analysis to have a file with snp id and position so let's write it in our folder 02_data
write.table(chr_pos, "02_data/SNP_pos.txt", sep="\t", quote=F, row.names=F)

#extract and format genotype matrix
geno <- extract.gt(obj.vcfR) # Character matrix containing the genotypes

#an empty matrix, (9 stands for missing data)
G <- matrix(9, nrow = nrow(geno), ncol = ncol(geno))

#that we fill with genotypes
G[geno %in% c("0/0", "0|0")] <- 0
G[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G[geno %in% c("1/1", "1|1")] <- 2

#an overview of our data and its first 10 rows/10 columns
table(as.vector(G))
dim(G)
G[1:10,1:10]

#As it will be useful later, I suggest that we export it now
write.table(G, "02_data/geno_matrix.txt", sep="\t", col.names=F, row.names=F)
```
We obtain a matrix of genotypes with 9 as missing data, 8018 rows for each SNP, and 240 columns for each individual
We will now ask outFLANK to calculate FST for eech locus. It needs the information about populations
For OutFLANK we will keep only the pop column. Then we will calculate a FST value for each SNP
```
#import pop info
info_samples_canada<- read.table("02_data/info_samples_canada.txt", header=T)
head(info_samples_canada)
pop_vector<- info_samples_canada$pop

# FST matrix with OutFLANK
my_fst <- MakeDiploidFSTMat(t(G), locusNames = id_snp, popNames = pop_vector)
```
Now we are ready to run OutFlANK.
We will follow the best practices recommended by the authors:
- remove SNPs with very low heterozygosity (options: Hmin = 0.1)
- use the FSt uncorrected for population size (options: NoCorr = TRUE) (anyway, here all pop have 20 individuals)
- Compare the FSt against a distribution based on independant SNPs (pruned for short-distance and long-distance LD)
We will use the list of pruned SNPs extracted with PLINK earlier.
Note that other possibilities exists such as using the package bigsnpr

####Step 1 we will run OutFLANK on the pruned SNPs and look at the distribution of FST
We will use the prune.in file produced by plink and a few manipulation to know the position of the SNPs that are in the pruned subset and will be include 

```
#import pruned info
id_snp_pruned<-read.table("02_data/population.can.random.snp.prune.in")

#those are SNPs id, we need to know at which position they are
#this can be done with the %in% function
lines_trim<-which(id_snp %in% id_snp_pruned[,1]) 
head(lines_trim)
length(lines_trim)

#run outFLANK on pruned SNPs
#numberOfSamples is the number of populations
#qthreshold is the false discovery rate
out_trim <- OutFLANK(my_fst[which(id_snp %in% id_snp_pruned[,1]),], NumberOfSamples=20, qthreshold = 0.05, Hmin = 0.1)
str(out_trim)

#have a look at the results
#the jpeg line allow to output an image in your folder that you can later download to have a look at
jpeg("04_outflank/outflank_prunedSNP_fst.jpeg")
OutFLANKResultsPlotter(out_trim, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =FALSE, RightZoomFraction = 0.05, titletext = NULL)
dev.off()

jpeg("04_outflank/outflank_prunedSNP_pvalues.jpeg")
hist(out_trim$results$pvaluesRightTail)
dev.off()
```

![outflank_trim1](07_img_readme/outflank_prunedSNP_fst.jpeg)
![outflank_trim1](07_img_readme/outflank_prunedSNP_pvalues.jpeg)
The p-value should be more or less flat and the distribution of FST about normal

####Step 2 We now run OutFLANK on all SNPs, corrected by the trim dataset
```
P1 <- pOutlierFinderChiSqNoCorr(my_fst, Fstbar = out_trim$FSTNoCorrbar, dfInferred = out_trim$dfInferred, qthreshold = 0.05, Hmin=0.1)
head(P1)

#we need to add the chromosome/position info for plotting. the left_join function in dplyr is super useful to match differebnt table
library(dplyr)
P1_pos<-left_join(P1, chr_pos, by=c("LocusName"="id_snp"))

#We can have a look at the results by exporting the figures
#we can look at the FSt as a function of heterozygosity to understand which snps have been evaluated, which one appear true or false outliers
And we can look along the genome with our manhattan plots

jpeg("04_outflank/outflank_outlier_fst_He.jpeg")
ggplot(P1_pos, aes(x=He, y=FST, colour=OutlierFlag))+ 
  geom_point()+
  theme_classic()
 dev.off()

#note that we divide here position by 1000000 so the scale is in MB
jpeg("04_outflank/outflank_outlier_fst.jpeg")
ggplot(P1_pos, aes(x=position/1000000, y=FST, colour=OutlierFlag))+ 
  geom_point()+
  theme_classic()+
  facet_grid(cols = vars(chromosome), scales = "free_x", space="free_x")+
  labs(  x = "position (in MB)")
 dev.off()

#It may also be easier to export the matrix and play in Rstudio
write.table(P1_pos, "04_outflank/outflank_fst_outliers.txt", sep="\t", row.names=F, quote=F)
```
![outflank_he](07_img_readme/outflank_outlier_fst_He.jpeg)
![outflank_along genome](07_img_readme/outflank_outlier_fst.jpeg)

What do you see? 
Lots of outliers are in chr5, and/or chr4... Are we really suprised?
Let's remember here that we are in a species with extremely high gene flow.  Moreover, we are looking here at outliers FST accross all populations (without giving any geographic or envrionmental information)
There is high differentiation between the two sexes and we have unbalanced sex-ratio in the sampling: this is probably driving the signal on chr 5. But there are also possibly adaptive loci on sex chromosomes... 

What about chr 4? It may be that the rearrangement is involved in local adaptation, or it is particularly divergent, and small fluctuations of frequency between populations are driving the signal?

Should we put a different threshold for FST outliers between different regions of the genome? And how to do it? This is an open question at the moment.
There is accumulating litterature suggesting that recombination is a very important factor, which may both create spurious effects, leading even under neutrality to high FSt variance but also, biological effets by promoting cluster of adaptive loci, without mentionning the matter of linked selection or recombination...

You may want to re-run on the SNP subset without chr4/chr5 to get a sense of what that says.

### With Baypass
To look at adaptive differentiation and environmental asssociations, we wil use Baypass http://www1.montpellier.inra.fr/CBGP/software/baypass/
The publication is here https://www.genetics.org/content/201/4/1555
Mathieu Gautier
GENETICS December 1, 2015 vol. 201 no. 4 1555-1579; https://doi.org/10.1534/genetics.115.181453
And there is a good manual here : http://www1.montpellier.inra.fr/CBGP/software/baypass/files/BayPass_manual_2.2.pdf

This package is an extension of Bayenv. It may be a little long to run so you will have time to skim the manual for more details. There are different underlying bayesian modle, and we will just run the default mode but you may want to explore a bit further in the future


#### Step 0: prepare files
We will use again scripts prepared by Yann Dorant in his toolbox:
```
python Toolbox/reshaper_baypass.py 02_data/population.can.random.snp.vcf 02_data/popmap_canada.txt 05_baypass/population.can.random.baypass
python Toolbox/reshaper_baypass.py 02_data/population.can.random.snp.pruned.recode.vcf 02_data/popmap_canada.txt 05_baypass/population.can.random.pruned.baypass
```

#### Step 1: run baypass basic

```
PATH_BP="your_path_to_the_source_of_baypass"
$PATH_BP/sources/g_baypass -npop 12 -gfile 05_baypass/population.can.random.baypass \
-outprefix 05_baypass/allsnps.output -nthreads 3
```
We can have a look at the output. the most importnat are the omega matrix, which is covariance between populations and the xtx value which represent a kind of value associated to each SNP. the higher it gets, the more this snps differentiate the populations.
Please copy those files so that we can look at them in Rstudio on your local computer while the next step will run on the server

#### Step 2: Have a threshold for the xtx value
While this runs, we can prepare our subsequent analysis. We need to know above which xtx threhold value we can consider a locus to be an outlier of population differentiation
To look for that, the authors suggest to simulate a neutral distribution with a small R function, run baypass on the simulated genotypes, and extract the distribution of XtX values. We can then chose a threshold of the 95% quantile, 99% quantile, etc... 

Open R in your terminal and let's simulate genotypes
```
args = commandArgs(trailingOnly=TRUE)

#source functions from Gautier file
source("$PATH_BP/utils/baypass_utils.R")

#load omega matrix
omega = as.matrix(read.table("05_baypass/prunedsnps.output_mat_omega.out"))

#load beta parameters
pi.beta.coef=read.table("05_baypass/allsnps.output_summary_beta_params.out", h=T)$Mean

#create the POD - simulate - don't forget to change suffix name
simu.bta <- simulate.baypass(omega.mat=omega, nsnp=2000,beta.pi=pi.beta.coef, suffix="simulates_pods")
```
You can now quit R, mv the simulated file into the folder baypass and run baypass on the simulated data 
```
mv G.simulates_pods 05_baypass/
$PATH_BP/sources/g_baypass -npop 12 -gfile 05_baypass/G.simulates_pods  -outprefix 05_baypass/simulate.output -nthreads 2
```
Please copy the output on the simulated dataset on your local computer.

#### Step 3: run baypass controlled for population structure
We want again to pinpoint the outlier that differ from the basic structure of the population. This is not a big deal in our example but it may be in other system. 
Best practices suggest to first run baypass on the LD-pruned vcf and to extract the covariance matrix between population (mat_omega), and then use it as a covariate in the baypass model

```
#run on pruned snps
$PATH_BP/sources/g_baypass -npop 12 -gfile 05_baypass/population.can.random.pruned.baypass \
-outprefix 05_baypass/prunedsnps.output -nthreads 2

#controlled run
"$PATH_BP"/sources/g_baypass -npop 12 -gfile 05_baypass/population.can.random.baypass \
-omegafile 05_baypass/prunedsnps.output_mat_omega.out \
-outprefix 05_baypass/allsnps.controlled.output -nthreads 2

#controlled run on simu dataset
$PATH_BP/sources/g_baypass -npop 12 -gfile 05_baypass/G.simulates_pods \
-omegafile 05_baypass/prunedsnps.output_mat_omega.out \
-outprefix 05_baypass/simulate_controlled.output -nthreads 2
```
Please copy the controlled output on your local computer.

#### A note
For best results, it is recommended to make several independant runs and take the median value for xtx or BF.
We dont have time and resource to do that today but you can try for your data with a loop and random starting seed
```
# run model for 5 repetitions
for i in 1 2 3 4 5 
do
seed=$((1000 + RANDOM % 9999)) #set random seed
echo "$seed"
g_baypass -npop $N_pop -gfile $geno -efile $env -omegafile $omega_mat -outprefix output"$i" -nthreads $N_CPU -npilot 25 -burnin 5000 -seed $seed
done
```

#### Visualising results on Rstudio on your local computer
If you are curious, you can explore the different file, look at the omega matrix, etc, which will be more or less like the pairwise fst matrix that we look at earlier
We will focus on xtx values
```
#load xtx values
xtx_allsnps<-read.table("05_baypass/allsnps.output_summary_pi_xtx.out", header = T)
head(xtx_allsnps)
#we will mostly work with M_WtX

#load position info about the SNPs.
SNP_pos<-read.table("02_data/SNP_pos.txt", header=T)
#should be same nb of rows
dim(xtx_allsnps)
dim(SNP_pos)

xtx_pos<-cbind(SNP_pos, xtx_allsnps)

ggplot(xtx_pos, aes(x=position, y=M_XtX, colour=chromosome))+ 
  geom_point()+
  theme_classic()+
  facet_grid(cols = vars(chromosome), scales = "free_x", space="free_x")
```
Now we realised that we really need to know at which value we put the threshold
```
#load xtx values from simulatd data
xtx_simu<-read.table("05_baypass/simulate.output_summary_pi_xtx.out", header=T)
head(xtx_simu)

#calculate the threshold
threshold_fdr0.01 = quantile(xtx_simu$M_XtX,probs=0.99)
threshold_fdr0.05 = quantile(xtx_simu$M_XtX,probs=0.95)

#add it on the plot
ggplot(xtx_pos, aes(x=position, y=M_XtX, colour=chromosome))+ 
  geom_point()+
  theme_classic()+
  facet_grid(cols = vars(chromosome), scales = "free_x", space="free_x")+
  geom_hline(aes(yintercept =threshold_fdr0.05), linetype="dotted", size=1, col="red", show.legend = FALSE)+
  geom_hline(aes(yintercept =threshold_fdr0.01), linetype="dotted", size=1, show.legend = FALSE)

#output outliers
xtx_pos[xtx_pos$M_XtX>= threshold_fdr0.05,]
```

Now you may want to do the same for the analysis controlled by population structure. Compare the results. What do you think? 
![xtx_uncontrolled](https://github.com/clairemerot/physalia_adaptation_course/blob/master/03_day3/07_img_readme/xtx_uncontrolled.png)
![xtx_controlled](https://github.com/clairemerot/physalia_adaptation_course/blob/master/03_day3/07_img_readme/xtx_controlled.png)

## 3C Environmental association
### With Baypass (locus by locus)
Baypass can also integrate environmental matrix or phenotypic matrix to run associations.
This means that we will no longer look at overall average differentiation between populations, but really ask which SNPs co-vary in frequency with environmental variables, in other words makes sens of differentiation between population along ecological gradients or look for phenotype/genotype association.

It is exactly the same procedure as before, except that we add the -efile option

#### step 0 format env file 
Baypass wants the environment as one line per env variable and one column per population. no header.
We can use our file 02_data/info_pop_geo_eco.txt to extract that. We will be interested in the temperature variable today.
As I did not pay attention to order the population alphabetically at the beginning of STACKS (sorry), we need to make really sure that our env file is exactly in the same order as the population are in the geno file. This is why there is a "order" column...

So please just go in R and try to output the env file well formatted for baypass.
R
```
#read file
info_pop<-read.table("02_data/info_pop_geo_eco.txt", header =T)
head(info_pop)
#env vector and its transpose form
info_pop$temperature
t(info_pop$temperature)
write.table (t(info_pop$temperature), "05_baypass/env.txt", sep="\t", quote=F, row.names=F, col.names=F)
quit()
```
Now we can run baypass with the environment file on both the real dataset and the simulated dataset.
Here is the example with and without the covariate matrix controlling for a possible underlying structure but depending on your system you may want to control or not, or do both to compare
```
$PATH_BP/sources/g_baypass -npop 12 -gfile 05_baypass/G.simulates_pods \
-efile 05_baypass/env.txt  \
-outprefix 05_baypass/simulate_env.output -scalecov -nthreads 3

"$PATH_BP"/sources/g_baypass -npop 12 -gfile 05_baypass/population.can.random.baypass \
-efile 05_baypass/env.txt  \
-outprefix 05_baypass/allsnps_env.output -scalecov -nthreads 2

$PATH_BP/sources/g_baypass -npop 12 -gfile 05_baypass/G.simulates_pods \
-efile 05_baypass/env.txt \
-omegafile 05_baypass/prunedsnps.output_mat_omega.out \
-outprefix 05_baypass/simulate_env.controlled.output -scalecov -nthreads 3

"$PATH_BP"/sources/g_baypass -npop 12 -gfile 05_baypass/population.can.random.baypass \
-efile 05_baypass/env.txt \
-omegafile 05_baypass/prunedsnps.output_mat_omega.out \
-outprefix 05_baypass/allsnps.controlled_env.output -scalecov -nthreads 2
```
#### Visualising results on Rstudio on your local computer
Like we did before we can plot the xtx but we will mostly be interested in the BF value (Bayesian factor of association with an environmental variable)
You will find it in the betai_reg.out in the column BF.dB.
"the Bayes Factor (column BF(dB)) in dB units (i.e., 10 × log10(BF)) measuring the support of the association of each SNP with each population covariable and the corresponding regression coefficients βi (column Beta_is)"
BF can be informative in itself, good candidate are usually above 20. We can also d find the threshold from the run on simulated data.
We will let you have a look at those data and plot BF along the genome re-using the same code as before.

What do you observe? 
Uncontrolled model
![bf_uncontrolled](https://github.com/clairemerot/physalia_adaptation_course/blob/master/03_day3/07_img_readme/bf_uncontrolled.png)

Controlled model
![bf_controlled](https://github.com/clairemerot/physalia_adaptation_course/blob/master/03_day3/07_img_readme/bf_controlled.png)

So on the data we have today, we find very few putative outlier SNPs associated with temperature. It may be that we don't have the power to detect some, we have only 12 populations and a reduced representation of the genome. It may be that given the high gene flow, very few regions of the genome are under strong enough selection to show detectable differences in allelic frequencies.

### With redundancy analysis (multi-loci)
So until now we have look at signatrue of adaptation locus by locus. This is a design that is more likely to point strong selective sweep, large-effect locus and large islands of divergence. Yet, most adaptation is expected to be polygenic and genome-scan are possibly not the best option to tackle this problem... Yet, doing RDa analysis which consider the whole genetic variance together with the geographic variation or ecological variation may be one of the first step in that direction.

Again, we can also used the RDA to identify putative candidate.
Today we will follow a tutorial inspired by this paper https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.14584 
Forester, BR, Lasky, JR, Wagner, HH, Urban, DL. Comparing methods for detecting multilocus adaptation with multivariate genotype–environment associations. Mol Ecol. 2018; 27: 2215– 2233. https://doi.org/10.1111/mec.14584

Which has a super nice vignette that you can check here https://popgen.nescent.org/2018-03-27_RDA_GEA.html

We will follow most steps with our dataset.
We will do most of the computational aspect of the RDA on the server but you could probably also run it locally, it is usually quite fast.

#### Step 0 Prepare data 
RDA takes a matrix of genotypes (0,1,2) with as many rows as samples and as many columns as loci. No missing data.We’ll use a simple approach to imputing missing values: we will impute using the most common genotype at each SNP across all individuals 

We can obtain following the same method as we did for the PCA yesterday. 
We can also use the geno matrix that we exported earlier from vcfR, and the associated SNP_pos file
```
geno<-read.table("02_data/geno_matrix.txt")
SNP_pos<-read.table("02_data/SNP_pos.txt", header=T)

#transpose data and give meaningful colnames
gen<-t(geno)
colnames(gen)<-paste(SNP_pos$chromosome, SNP_pos$position, sep="_")
gen [1:10,1:10]

#replace 9 by NA
gen[which(gen=="9")]<- NA
#evaluate % of missing
sum(is.na(gen))/(dim(gen)[1]*dim(gen)[2]) # about 1.5% of missing data
#impute missing with the most common geno
gen.imp <- apply(gen, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
sum(is.na(gen.imp)) # No NAs

```
Let's look now at our environmental/pheno matrix. 
```
info<-read.table("02_data/info_samples_canada.txt", header=T)
head(info)
```


#### Step 1 Run rda on environmental variable and test it 
The command lines to run a rda are quite simple:

#run rda
```
library(vegan)
temp.rda <- rda(gen.imp ~ info$temperature, scale=T)
temp.rda
```
Since we have only one variable, only RDA1 (the 1st axis is meaningful)
We can look at the fraction of variance explained by this RDA and its significance using permutation tests. 

```
RsquareAdj(temp.rda)

signif.full <- anova.cca(temp.rda, parallel=getOption("mc.cores")) # default is permutation=999
signif.full
```
so you see that temperature does not explain a large fraction of total variance (0.5%!!) but the model appear significant. It likely means that very few loci covary with temperature but this parameter does not explain a large fraction of genetic variation as a whole.


#### Step 2 Analyse RDA output







#### Step 3 environmental variation vs.geographic variation?
With the RDA we can also ask about the relative contribution of geography and environment since the rda can be controlled by a third matrix. It is also possibel to include a previously inferred geentic structure, etc.

Very simply we can use latitude and longitude

```
geo.rda <- rda(gen.imp ~ info$lat + info$longe, scale=T)
geo.rda
RsquareAdj(geo.rda)

temp.geo.rda <- rda(gen.imp ~ info$temp, cbind(info$lat,info$longe), scale=T)
temp.geo.rda
RsquareAdj(temp.geo.rda)

signif.axis <- anova.cca(geo.rda, by="axis", parallel=getOption("mc.cores"))
signif.axis
```

```
#Create a geographic matrix with dbmem
coorgeo<-info[,6:7]
dist_eucl<-daisy(coorgeo)
geo_pcnm<-pcnm(dist_eucl)
geo_mat<-geo_pcnm$vectors
```

#### Repeat the analysis on phenotype?


#### expand to multivariate?
Another powerful aspect of the RDA is to be able to analyse together several environmental variables. 

Yet, they should not be correlated. so if you include several environmental axis, it is best to first do quick correlation steps and select the most interesting variables or PC axis that summarize them. 

The problème of including several variables (or worse sort of composite variables based on PCA is that it can become tricky to disentangle which outliers are associated with which environmental variable). 

Today, we only included temperature but if you wish, you can use the gps points and databases (bioclim, marspec...) to retrieve more environmental data. If you are interested into that we can come back on it later or on the last day. We palso resent here a sort of dummy example with latitude and longitude as variables if you have time to run the multivariate.

If you include several environmental variables (for instance rda(gen.imp~info$temperature +info$longitude), you can test of teh association between genetic variation and each variable with the model by axis, and check that the correlation between environmental variable with the VIF. It might take a bit longer to run 

```
temp_gps.rda <- rda(gen.imp ~ info$temperature+ info$latitude + info$longitude, scale=T)
temp_gps.rda

RsquareAdj(temp_gps.rda)

signif.axis <- anova.cca(temp_gps.rda, by="axis", parallel=getOption("mc.cores"))
signif.axis
vif.cca(temp_gps.rda)

```

