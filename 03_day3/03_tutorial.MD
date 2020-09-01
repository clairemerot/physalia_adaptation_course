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
We will calculate here FST between all pairs of populations. FST accross the genome are expected to be largely driven by neutral markers while peaks of FSt may be related to selection for local adaptation... So we will first try to get a sense of global pattern and then look for outliers

to get pairwise FSt, we will use STAMP with a toolbox developped by Yann
```
```
Now you can export the pairwise FST matrix on your local computer to make a few visualisations in R
Please export the info files about populations too

We are now in Rstudio on your computer.
We keep it simple and do a simple numeric matrix but you can imagine more fancy ways, with heatmaps or so.

```
library(corrplot)
corrplot?
```

What do you notice? Is it heterogeneous? Do some population look more differentiated than other?

### Isolation by distance
To explore whether this is linked to distance between population, we will do a IBD test (Isolation by distance)
```
library(SoDA)
#import information about populations
info_pop<-read.table("02_data/info_pop_geo_eco.txt", header=T)
head(info_pop)
#calculate geogrpahic (euclidian) distances between all pairs of populations
distance<-dist(geoXY(info_pop$latitude, info_pop$longitude))
distance_matrix<-as.matrix(distance) #make a triangular matrix
colnames(distance_matrix)<-gps$pop
rownames(distance_matrix)<-gps$pop
distance_matrix

#make a tabular matrix
dist_tab<-vector(length=0)
for (i in 1:12)
{
  dist_tab<-c(dist_tab, distance_matrix[i, (i+1):13])
}

#test association with FST

#plot IBD

```

does it make sense? 
-> I hope it does not. 
If populations that have only males looks closer we can make them think about sex-linked markers and the importance of sampling bias. 
Can we look at the FSt landscape along the genome for one pair of pop? Let's say for the most differentiated pop?If peaks in chr5 -> need to pay attention to sex


Let's look at the sex ratio in our data.
```
info_ind<-read.table("02_data/info_samples_canada.txt", header=T)
head(info_ind)
table (info_ind$pop,info_ind$sex)
```
It seems that the field sampling has not been very good at balancing sex-ratio between populations... We should be worried about sex-linked markers driving the signal!
(note: we intentionnaly subsetted the dataset to create this bias which was not in the original publication ;-) but it may happen easily for some species or low sample size)

### Getting rid of LD spurious effects?
Yesterday, we saw that a region of chromosome 4 and sex-linked markers on chr 5 were overwhelming the structure. Will that influence our pairwise estimates? Possibly
Let's re-run the steps above on the vcf in which we removed the chr5 and half of the chr4.
Let's look at the pairwise FST matrix and the IBD stats
What do you see now?

### Making a LD-pruned vcf
Another way to get rid of SNPs in LD due to the inversion or sex (?) would be to use a LD-pruned set of SNPs. For this we can use PLINK which will compute LD between SNPs by windows along the genome and keep one SNP out of several SNPs in linkage. Here we are interested in long distance LD (since we have already removed short-distance LD by keepin only one SNP per RAD locus) so I have set the window quite large.

Claire add plink code
```
plink
```

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
obj.vcfR <- read.vcfR("02-data/population.can.random.snp.vcf")

#following the tutorial
geno <- extract.gt(obj.vcfR) # Character matrix containing the genotypes
position <- getPOS(obj.vcfR) # Positions in bp
chromosome <- getCHROM(obj.vcfR) # Chromosome information
chr_pos<-paste(chromosome, position, sep="_") # make unique locus name

#an empty matrix, (9 stands for missing data)
G <- matrix(9, nrow = nrow(geno), ncol = ncol(geno))

#that we fill with genotypes
G[geno %in% c("0/0", "0|0")] <- 0
G[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G[geno %in% c("1/1", "1|1")] <- 2

#an overview of our data
table(as.vector(G))
dim(G)
```
We obtain a matrix of genotypes with 9 as missing data, 8018 rows for each SNP, and 240 columns for each individual
We will now ask outFLANK to calculate FST for eech locus. It needs the information about populations
For OutFLANK we will keep only the pop column. Then we will calculate a FST value for each SNP
```
#import pop info
info_samples_canada<- read.table("02-data/info_samples_canada.txt", header=T)
head(info_samples_canada)
pop_vector<- info_samples_canada$pop

# FST matrix with OutFLANK
my_fst <- MakeDiploidFSTMat(t(G), locusNames = chr_pos, popNames = pop_vector)
```
Now we are ready to run OutFlANK.
We will follow the best practices recommended by the authors:
- remove SNPs with very low heterozygosity (options: Hmin = 0.1)
- use the FSt uncorrected for population size (options: NoCorr = TRUE) (anyway, here all pop have 20 individuals)
- Compare the FSt against a distribution based on independant SNPs (pruned for short-distance and long-distance LD)
We will use the list of pruned SNPs extracted with PLINK earlier.
Note that other possibilities exists such as using the package bigsnpr

####Step 1 we will run OutFLANK on the pruned SNPs and look at the distribution of FST
```
#import pruned info
which_pruned<-read.table()

#run outFLANK on pruned SNPs
#numberOfSamples is the number of populations
#qthreshold is the false discovery rate
out_trim <- OutFLANK(my_fst[which_pruned,], NumberOfSamples=20, qthreshold = 0.05, Hmin = 0.1)
str(out_trim)

#have a look at the results
#the jpeg line allow to output an image in your folder that you can later download to have a look at
jpeg("04-outflank/outflank_prunedSNP_fst.jpeg")
OutFLANKResultsPlotter(out_trim, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =FALSE, RightZoomFraction = 0.05, titletext = NULL)
dev.off()

jpeg("04-outflank/outflank_prunedSNP_pvalues.jpeg")
hist(out_trim$results$pvaluesRightTail)
dev.off()
```
The p-value should be more or less flat and the distribution of FST about normal

####Step 2 We now run OutFLANK on all SNPs, corrected by the trim dataset
```
#run OutFLANK
P1 <- pOutlierFinderChiSqNoCorr(my_fst, Fstbar = out_trim$FSTNoCorrbar, dfInferred = out_trim$dfInferred, qthreshold = 0.05, Hmin=0.1)
head(P1)

#add the chromosome/position info


#We can have a look at the results by exporting the figures
jpeg("04-outflank/outflank_outlier_fst.jpeg")
ggplot(P1, aes(x=position, y=FST, colour=OutlierFlag))+ 
  geom_point()+
  theme_classic()+
  facet_grid(cols = vars(chrom), scales = "free_x", space="free_x")
 dev.off()

#Or export the matrix and play in Rstudio
write.table(P1, "04-outflank/outflank_fst_outliers.txt", sep="\t", row.names=F, quote=F)
```


### With Baypass


## 3C Environmental association
### With Baypass (locus by locus)

### With redundancy analysis (multi-loci)
