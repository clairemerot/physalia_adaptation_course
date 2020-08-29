# Day2: Population structure and differentiation
## Step 1: Calculating F-statistics
### With Stacks
As you can see if you ```cd``` in one of the folders storing the output from your ```populations``` runs, STACKS produces many output files. The manual provides a thorough explanation of all the different file formats. Go to <https://catchenlab.life.illinois.edu/stacks/manual/> and scroll down to section 6.6. For each output type, open the corresponding file on the server, and take your time to explore them.

One statistics we are interested in F<sub>ST</sub>, a relative measure of differentiation between populations. Let's look at one of these files in more detail
```
cd ~/stacks/populations_2lin_random
head -n 30 populations.fst_NWA-GRE.tsv
```
F<sub>ST</sub> for each locus is at the 8<sup>th</sup> column. We can look at the 20 most differentiated SNPs with
```
grep -v "^#" populations.fst_NWA-GRE.tsv | sort -rk8,8 | head -n 20
```
```grep -v``` exclude lines that start with #, then you sort the file by values of F<sub>ST</sub> from the highest to the lowest (```-r```), and print to the screen the top 20 SNPs.

If you have more than two populations, STACKS will print out ```populations.fst_POP1-POP2.tsv``` for each unique combination of populations. Look into the ```populations_canada_random``` folder. Here we can create a list of the 100 most differentiated SNPs among all populations with this one-liner
```
cat populations.fst_*.tsv | grep -v "^#" - | sort -rk8,8 | head -n 1000 | sort -un -k1,1 -k2,2 -s | sort -rk8,8 | cut -f 1,6 | tail -n +2 | head -n 100 | sort -n > high_fst.whitelist.tsv
```
which is mostly based on commands that we've just used above, with the difference that you merge all files containing population pairwise information and print only the SNP information. Note that even though the data are mapped to a reference genome, STACKS has built a catalog of loci and uses locus ID and position within the locus (stored in columns 1 and 6 in this file) to retrieve information on these SNPs (i.e. STACKS doesn't use chromosome and position to identify SNPs).

Have a look at this list. Are these SNPs randomly distributed across the 5 chromosomes?

Now that we have created a 'whitelist', we can rerun populations and generate F-statistics and input files for these 100 highly differentiated SNPs only.
```
cd ~/stacks
mkdir populations_canada_random_highfst
populations -t 16 -P ~/stacks/gstacks/ -M ~/scripts/popmap_canada.txt -O populations_canada_random_highfst --fstats --vcf --genepop --structure -W ~/stacks/populations_canada_random/high_fst.whitelist.tsv
```


Whitelists and blacklists are very useful to select or discard loci for analysis. For quick data exploration, or for analysis that don't require, or can't handle, large datasets, you can use a similar code to create a whitelist with a random subset of SNPs. For 1000 SNPs for example
```
grep -v "^#" populations.sumstats.tsv | cut -f 1,4 | sort | uniq | shuf | head -n 1000 | sort -n > 1000snps_whitelist.tsv
```
Note that STACKS keep into account the whole locus, rather than just the random SNP we selected in each locus, for many statistics. Open ```populations.sumstats_summary.tsv``` for example. To help visualization you can download the file and import it in excel. Additionally, STACKS provides many haplotype-based statistics. Although we won't delve into these today, these are important to obtain estimates of D<sub>xy</sub>, a measure of absolute divergence, gene and haplotype diversity, and to phase SNPs, which can provide additional important information for both evolutionary questions and conservation applications.

See following review on the use of haplotype information in conservation genomics
```
Leitwein, M., Duranton, M., Rougemont, Q., Gagnaire, P.A. and Bernatchez, L., 2020.
Using haplotype information for conservation genomics.
Trends in Ecology & Evolution, 35(3), pp.245-258.
```

Note also that statistics will be different depending on whether you look at all SNPs in one locus or only one random SNP. See below

Statistics from ~/stacks/populations_2lin/populations.log, which includes all SNPs
```
Population summary statistics (more detail in populations.sumstats_summary.tsv):
  NWA: 39.172 samples per locus; pi: 0.054733; all/variant/polymorphic sites: 770730/45238/25848; private alleles: 18011
  GRE: 39.138 samples per locus; pi: 0.054095; all/variant/polymorphic sites: 770730/45238/27188; private alleles: 19351

Population pair divergence statistics (more in populations.fst_summary.tsv and populations.phistats_summary.tsv):
  NWA-GRE: mean Fst: 0.024875; mean Phi_st: 0.077792; mean Fst': 0.074473; mean Dxy: 0.0049455
```
Statistics from ~/stacks/populations_2lin_random/populations.log, which includes only one random SNP/locus.
```
Population summary statistics (more detail in populations.sumstats_summary.tsv):
  NWA: 39.161 samples per locus; pi: 0.051743; all/variant/polymorphic sites: 770730/8530/4634; private alleles: 3345
  GRE: 39.12 samples per locus; pi: 0.054156; all/variant/polymorphic sites: 770730/8530/5171; private alleles: 3882

Population pair divergence statistics (more in populations.fst_summary.tsv and populations.phistats_summary.tsv):
  NWA-GRE: mean Fst: 0.028809; mean Phi_st: 0.033052; mean Fst': 0.020139; mean Dxy: 0.00092926
```

Here, the value of π are not very different because they are averaged over only the polymorphic sites (which is NOT how we should estimate π anyway). Same for F<sub>ST</sub>. Do you see how different the haplotype-based estimates are? They are all lower in the random-SNP dataset because STACKS computes these statistics as if there was only one polymorphic site in each of these loci.


### With VCFtools
To calculate F<sub>ST</sub> in VCFtools you need to provide the list of samples for each population to include in the calculation. For the reduced dataset, including only 80 individuals from Greenland and Canada, you will find the list of samples for each of the two lineages in the ```~/scripts``` folder
```
vcftools --vcf populations_2lin_random/populations.snps.vcf --weir-fst-pop ~/scripts/pop_canada.txt --weir-fst-pop ~/scripts/pop_greenland.txt  --out fst_2lin
```
As stdout you get
```
VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf populations_2lin_random/populations.snps.vcf
	--weir-fst-pop /home/ubuntu/scripts/pop_canada.txt
	--weir-fst-pop /home/ubuntu/scripts/pop_greenland.txt
	--keep /home/ubuntu/scripts/pop_canada.txt
	--keep /home/ubuntu/scripts/pop_greenland.txt
	--out fst_2lin

Keeping individuals in 'keep' list
After filtering, kept 80 out of 80 Individuals
Outputting Weir and Cockerham Fst estimates.
Weir and Cockerham mean Fst estimate: 0.032662
Weir and Cockerham weighted Fst estimate: 0.23741
After filtering, kept 8530 out of a possible 8530 Sites
Run Time = 0.00 seconds
```
The calculation is very fast, and in addition to per-site estimates you get mean and weighted F<sub>ST</sub> estimates, which you could save in a separate file by adding ```> fst_2lin.log``` at the end of the previous command.

If you have a reference genome and panel of SNPs dense enough (or whole genome resequencing data) you can calculate window estimates with
```
vcftools --vcf populations_2lin_random/populations.snps.filter.vcf --weir-fst-pop ~/scripts/pop_canada.txt --weir-fst-pop ~/scripts/pop_greenland.txt  --out fst_2lin_win --fst-window-size 100000 --fst-window-step 100000
```
In fact, window estimates are a practical way to summarize data and to reduce the number of datapoints to handle, which can be problematic for visualization.

We can also calculate π in VCFtools with the command
```
cd ~/stacks
vcftools --vcf populations_2lin_random/populations.snps.vcf --site-pi --keep ~/scripts/pop_canada.txt --out pi_canada
```
However, these single-site estimates of π are not very meaningful. π is defined as the average number of nucleotide differences per site between two DNA sequences in all possible pairs in the sample population. That is the same reason why π averaged across polymorphic sites in STACKS is not meaningful either. π must be averaged across the whole length of the sequence you are analyzing, whether it is a collection of RAD loci or whole chromosomes.

#### Important! There are many different ways to estimate FST and you should always report the program and the methods used.
For example, STACKS adopts an AMOVA F<sub>ST</sub> by Weir (from Genetic Data Analysis II, chapter 5, "F Statistics," pp166-167) and VCFtools uses the F<sub>ST</sub> estimate by Weir and Cocherham (1984). These two, as you can see from the analyses above, will give you slightly different estimates. Furthermore, VCFtools provided a weighted estimate of F<sub>ST</sub> that may be better suited to capture genome-wide population differentiatiation. 

## Step 2. Assessing population structure
### Structure with faststructure
STRUCTURE has been the most popular software for detecting and describing population structure since it was published 20 years ago. However, it was developed before we started genotyping thousands of markers at a time and it is therefore slow for large genomic datasets. fastSTRUCTURE (Raj eet al. 2014, Genetics) provides a much faster alternative for these datasets while it is based on a similar Bayesian framework.

We'll use the ```.structure``` file created by ```populations```. The first 6 columns of the file will be ignored; these typically would include IDs, metadata, etc. so we need to have 6 dummy columns before the genotypes for each individual/locus. Note that this software only handles bi-allelic loci. The two alleles at each locus can be encoded as desired; however, missing data should be encoded as '-9' (also check https://rajanil.github.io/fastStructure/ for details). Below is a little script to prepare the file and run fastSTRUCTURE.
```
cd  populations_all_random ### or populations_2lin_random or populations_canada_random
tail -n +3 populations.snps.structure > populations.structure_nohead ### deletes first two rows
cut -f1-2 populations.structure_nohead > twocolumns
paste twocolumns twocolumns populations.structure_nohead | sed 's/\t0/\t-9/g' > populations.faststructure.str
python structure.py -K 2 --input=populations.faststructure --output=capelin.all --format=str
```
Now you can modify the ```-K``` parameters and test different different numbers of clusters.
With this script you can test what number of populations/cluster best describes your dataset.
```
python chooseK.py --input=capelin.all
```
Finally, we can plot the admixture proportions with
```
python distruct.py

Here is how you can use this script

Usage: python distruct.py
     -K <int>  (number of populations)
     --input=<file>  (/path/to/input/file; same as output flag passed to structure.py)
     --output=<file>   (/path/to/output/file)
     --popfile=<file>  (file with known categorical labels; optional)
     --title=<figure title>  (a title for the figure; optional)

python distruct.py -K 2 --input=capelin.all --output=capelin.all_distruct
```
### PCA (Yann)

**[YANN] blabla about PCA...I'm writitng some few lines in draft...**



For this analysis, we will convert vcf files into geno format (matrix of genotypes) using the programm **vcftools**.
Here, we will work with two vcf files
* vcf file containing four populations (H,L,O and U)
* vcf file containing 12 populations (A-L)

So within your own folder in the Amazon server, use this two command (cmd) lines:
```
vcftools --vcf stacks/populations_2lin_random/populations.snps.vcf --012 --out populations_2lin_random
vcftools --vcf stacks/populations_can_random/populations.snps.vcf --012 --out populations_can_random

```
One this is done, you will see via the cmd line `ls` (for listing all files in the directory), that six new files have been generated.
* files ending with only `.012` contain genotypess information (0 for Ref homo, 1 for het, 2 for Alt homo and -1 for misisng data)
* files ending with `.012.pos` contain a list of SNP ids (Chromosme and Position) from the vcf file
* files ending with `.012.indv` contain a list of samples ids from the vcf file

So to work with this files, we will return in your local computer and download them.
Before downloading these files in your local computer, open a new terminal and move into your current working folder for this course ``PATH/physalia/``.
Then, you can download these files via the following command line :
```
scp -i ./boh.pem username@xx.xx.xx.xx:./*.012* ./02_day2/02-data/
```

## Here I think we need to shape a unique folder architecture because it could be annoying to deal with multiples files which are misstored.

Once the six files have been downloaded, open R studio and set your working directory to ``PATH/physalia/02_day2/``

#### Step 0: load the required libraries.
```
#Libraries
  library(dplyr)
  library(magrittr)
  library(tibble)
  library(ggplot2)
```
#### Step 1: load initial required data files
```
#1. load population map
popmap <- read.table("../00_documents/info_samples.csv", h=T,sep=';')

#2. load geno data for the 2lin
geno.012_2lin.012 <- read.table('02-data/population.2lin.rand.snp.012')[,-1] #load genotype matrix
geno.012_2lin.012.pos <- read.table('02-data/population.2lin.rand.snp.012.pos') %>% #load SNPs info
  mutate(., locus=paste(V1,V2,sep='_')) #create a new column for SNP info name (CHR + position)
geno.012_2lin.012.indv <- read.table('02-data/population.2lin.rand.snp.012.indv') #load individuals info

Set rownames and colnames to the geno matrix
dimnames(geno.012_2lin.012) <- list(geno.012_2lin.012.indv$V1, geno.012_2lin.012.pos$locus)
#check the geno matrix
geno.012_2lin[1:6,1:6]
```
```
Chr1_4669 Chr1_37123 Chr1_53559 Chr1_64218 Chr1_73377 Chr1_76193
L_01         0          0          0          0          0          1
L_02         0          0          0          0          0          0
L_03         0          0          0          0          0          0
L_04         0          0          0          0          0          0
L_05         0          0          0          0          0         -1
L_06         0          0          0          0          0          0
```
#### Step 2: Impute missing data
Missing data in 012.geno files from vcftools are coded with -1. We will change it for NAs
```
geno.012_2lin.012[geno.012_2lin.012 == -1] <- NA
```
```
Chr1_4669 Chr1_37123 Chr1_53559 Chr1_64218 Chr1_73377 Chr1_76193
L_01         0          0          0          0          0          1
L_02         0          0          0          0          0          0
L_03         0          0          0          0          0          0
L_04         0          0          0          0          0          0
L_05         0          0          0          0          0         NA
L_06         0          0          0          0          0          0
```

Here we will fill the NAs values by the most common genotype across all samples for a given SNP
```
geno.012_2lin.012.imp <- apply(geno.012_2lin.012,2,function(x){
                           replace(x, is.na(x), as.numeric(names(which.max(table(x))))) })
```
```
Chr1_4669 Chr1_37123 Chr1_53559 Chr1_64218 Chr1_73377 Chr1_76193
L_01         0          0          0          0          0          1
L_02         0          0          0          0          0          0
L_03         0          0          0          0          0          0
L_04         0          0          0          0          0          0
L_05         0          0          0          0          0          0
L_06         0          0          0          0          0          0
```

In genomics, PCA are multivariate analyses which could be bias by ultra rare variants (i.e. only represented by one sample)
>Linck, E., & Battey, C. J. (2019). Minor allele frequency thresholds strongly affect population structure inference with genomic data sets. Molecular Ecology Resources, 19(3), 639–647. doi: 10.1111/1755-0998.12995
https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.12995


Based on this avenue, we will filter our geno matrix by keeping only SNPs where the alternative allele is represented by at least two samples.
Usually, we name this filter as MAS for minimum number of samples with rare allele <int> eg: 2 or more.

Here I give you a custom R function to use this filter:
```
MAS <- function(geno, MAS.thresh =2){
#count number of samples which are heterozygous or rare homozygous
  get_MAS_locus <- apply(geno,2,function(x) sum(table(x)[-1]))
#get the list of locus that did not parse the threshold
  blacklist <- get_MAS_locus[get_MAS_locus < MAS.thresh]
#Print number of SNPs filtered out
  message(length(blacklist)," SNPs removed")
#return the filtered matrix
  return(geno[,-which(colnames(geno) %in% names(blacklist))])
}
```
To use it, this is simple.
```
geno.012_2lin.012.imp.MAS2 <- MAS(geno=geno.012_2lin.012.imp, MAS.thresh=2)
```
**Question : How many SNPs did not passed the MAS filter ?**

Next, we are ready to perform the PCA. The code is very simple. Only one line !
```
pca.2lin <- prcomp(geno.012_2lin.imp.MAS2, scale=T)
```
**(YANN) talk about PCA summary stats...**

OK ! now we are ready to make a plot of our PCA. To do this, there is so many ways and R libraries.
I will present you my own with the R library ggplot2, which is a highly modulable way.

First, we have to make a synthetic dataframe which incorporate various features:
* Number of PCs that we want to examine. Usually I keep the first four PCs, but you can keep more.
* Sample ids
* Popupations info
* any other caracteristic that you want related to the samples (e.g. size, color, sex...)

So how we do that ?
```
#prepare dataset to plot PCAs
pca.2lin.sub <- pca.2lin$x[,1:4] %>% #retain the first four PCs
  as.data.frame(.) %>% #transform to dataframe object
  tibble::rownames_to_column(., var='id') %>% #set rownames to column samples ids
  dplyr::left_join(., popmap, by='id') #Here we use the left_join function
                                       #from dplyr to wrap the population vector
                                       #of our samples.
```
Almost done! Now we will use the ggplot library to make a super figure !
```
ggplot(pca.2lin.sub) + aes(x=PC1, y=PC2, col=pop) +
  geom_hline(yintercept = 0, lty=2, col='grey50') + #add horiz line at y=0
  geom_vline(xintercept = 0, lty=2,col='grey50') +  #add vertical line at x=0
  geom_point() + #add the samples
  scale_color_manual(values=c('#cab2d6','purple','royalblue', 'cyan')) + #define a new color scale
  theme_bw() + #use  classic dark-on-light ggplot2 theme
  theme(panel.background = element_rect(fill='white'), #set some theme options
        panel.grid = element_blank())

#Print it !
ggsave("03-analyses/02-PCA/PCA_biplpot_2lin_capelin.png", width = 6, height = 5)
```
Done ! Good job.
![PCA_img](04-readme_img/PCA_1_tutorial_capelin.png)

Now, the question is what do you see ?

* **NOte:** Don't forget to save your script in your own local folder ``02_day2/00-scripts/``

Ok, that was easy no ? Well, I'm sure that you are able to do the same work for the other dataset composed by the 12 populations.
Test it and tell us what are your conclusions about it ?

I think you should be able to do this in <20 minutes.




### DAPC
A Discriminant Analysis of Principal Components is a multivariate approach that merges a Principal Component Analysis (PCA) and a Discriminant Analysis (DA). A PCA aims to summarize the variation among individuals and it runs very fast, also on large datasets. However, it is not powerful at discriminating groups because it doesn't use any a priori information on grouping, and intra-group variation can overwhelm inter-group variation. A DA, on the other hand, tries to summarize the variation among groups, while minimizing the variatin within groups. Thus, a DAPC takes the best of the two analyses to describe population structure. Also, compared to structure-like analysis it is not based on strict model assumptions and is more powerful at describing isolation-by-distance and hierarchical structure that the Bayesian approaches implemented in STRUCTRURE, fastSTRUCTURE, or ADMIXTURE.

On your desktop, create three new folders named ```populations_2lin_random```, ```populations_all_random``` and  ```populations_canada_random``` and download the file ```populations.snps.structure``` from each of three corresponding ```populations``` folders on the server. Change extension ```.structure```to ```.str``` of each of three files and open R Studio. Let's start analyzing the reduced dataset including 80 individuals from Canada and Greenland.
```
### load adegenet package
library(adegenet)
### set working directory
setwd("~/Desktop")
### load dataset and convert it from structure to genind format
twolin<-import2genind("populations_2lin_random/populations.str") # 80 ind and 8530 SNPs in my case, check your populations.log file to get the right number of SNPs
### This command will prompt questions about the structure of the file
```
![convert4dapc](https://github.com/clairemerot/physalia_adaptation_course/blob/master/images_tutorial/adegenet_loaddata.png)
```
### Once the dataset is loaded, you can proceed with the DAPC
dapc_twolin<-dapc(twolin) ### choose 50 PCs and 1 distriminant function, as you did have a choice...
###To plot results
scatter(dapc_twolin)
```
Because our dataset contains only two groups, NWA (Canada) and GRE (Greenland), we have only one discriminant function available. For this reason, we can plot our results across only one axis of variation. However, you can see that that's enough to sharply and unequivocally separate the Canadian and Greenlandic lineages.

You can repeat this analysis with the full dataset ```populations_all_random/populations.str```.
```
all<-import2genind("populations_all_random/populations.str") # 280 ind and 7943 SNPs
dapc_all<-dapc(all) ### play with number of PC and DF
scatter(dapc_all)
```
What do you see now?

If any structure is present among the Canadian populations, it may be hidden by the strong differentiation between the Canadian and Greenlandic lineages.
So let's repeat the analysis including only the Canadian populations, and use this to explore how to select the right number of PCs. Run and visualize the DAPCs based on 100 and 200 PCs.
```
canada<-import2genind("populations_canada_random/populations.str") # 240 ind and 8080 SNPs
dapc_canada1<-dapc(canada, n.da=4, n.pca=100)
dapc_canada2<-dapc(canada, n.da=4, n.pca=200)
scatter(dapc_canada1)
```
![dapc_canada1](https://github.com/clairemerot/physalia_adaptation_course/blob/master/images_tutorial/plot_dapc_canada_100_4.png)
```
scatter(dapc_canada2)
```
![dapc_canada2](https://github.com/clairemerot/physalia_adaptation_course/blob/master/images_tutorial/plot_dapc_canada_200_4.png)


If too few PCs (with respect to the number of individuals) are retained, useful information will be excluded from the
analysis, and the resultant model will not be informative enough to accurately discriminate
between groups. By contrast, if too many PCs are retained, this will have a destabilising
effect on the coefficients extimated, leading to problems of overfit. In such cases, the
model is able to describe all of the data in such detail that it becomes flexible enough
to discriminate almost perfectly between any possible clusters.

However, we can assess the trade-off between power of discrimination and over-fitting by calculating the alpha-score, which is the difference between the proportion of successful reassignment of the analysis (observed discrimination) and values obtained using random groups (random
discrimination).
```
temp1 <- optim.a.score(dapc_canada1)
temp2 <- optim.a.score(dapc_canada2)
```
These analyses suggest I should use fewer PCs, around 50.
```
dapc_canada1<-dapc(canada, n.da=4, n.pca=50)
```
![dapc_canada2](https://github.com/clairemerot/physalia_adaptation_course/blob/master/images_tutorial/plot_dapc_canada_50_4.png)


Another way to see the effect of the choice of the number of PCs is from stucture-like plots made with the function ```compoplot```, which plots the assignemnt proportions of each individual.
```
compoplot(dapc_canada1)
```
![compoplot1_canada1](https://github.com/clairemerot/physalia_adaptation_course/blob/master/images_tutorial/compoplot_dapc_canada1.png)
```
compoplot(dapc_canada2)
```
![compoplot1_canada1](https://github.com/clairemerot/physalia_adaptation_course/blob/master/images_tutorial/compoplot_dapc_canada2.png)
```
compoplot(dapc_canada3)
```
![compoplot1_canada1](https://github.com/clairemerot/physalia_adaptation_course/blob/master/images_tutorial/compoplot_dapc_canada3.png)


Although 200 PCs discriminate populations well, we know from low pairwise population F<sub>ST</sub> that differentiation is very low, which lead us to conclude that these high population assignments are the consequence of data over-fitting. On the other hand, the poor population assignemnt you obtain from 50 PCs is indicative of the weak population structure among the Canadian populations.

You can change pretty much everything in your plot and the authors of adegenet have put together a great tutorial <https://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf>

## Step 3. Assessing population structure heterogeneity across the genome
### Sliding PCA (Claire)
As you saw earlier, the PCA we obtained is very unexpected and it looks like some portion of the genome are overwhelmingly driving the structure.
To get a better sense of what's going on, we will be doing PCA again, but along the genome by window of X SNPs.
For this we will use a R package available here https://github.com/petrelharp/local_pca 
And presented in this publication https://www.genetics.org/content/211/1/289


### Repeat analyses on different genomic regions/chromosomes
### Manhattan plots for sex and inversion
