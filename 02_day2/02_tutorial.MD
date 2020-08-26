# Day2: Population structure and differentiation 
### Step 1: Calculating F-statistics
#### With Stacks
As you can see if you ```cd``` in one of the folders storing the output from your ```populations``` runs, STACKS produces many output files. The manual provides a thorough explanation of all the different file formats. Go to <https://catchenlab.life.illinois.edu/stacks/manual/> and scroll down to section 6.6. For each output type, open the corresponding file on the server, and take your time to explore them.

One statistics we are interested in F<sub>ST</sub>, a relative measure of differentiation between populations. Let's look at one of these files in more detail
```
cd ~/stacks/populations_canada_random
head -n 30 populations.fst_L-K.tsv
```
F<sub>ST</sub> for each locus is at the 8<sup>th</sup> column. We can look at the 20 most differentiated SNPs with
```
grep -v "^#" populations.fst_L-K.tsv | sort -rk8,8 | head -n 20
```
```grep -v``` exclude lines that start with #, then you sort the file by values of F<sub>ST</sub> from the highest to the lowest (```-r```), and print to the screen the top 20 SNPs.

We can also create a list of the 100 most differentiated SNPs among all populations with this one-liner
```
cat populations.fst_*.tsv | grep -v "^#" - | sort -rk8,8 | head -n 200 | cut -f 1,6 | sort | uniq | sort -n | head -n 100 > high_fst.whitelist.tsv
```
that is mosyly based on commands that we've just used above, with the difference that you merge all files containing population pairwise information and print only the SNP information. Note that even though the data are mapped to a reference genome, STACKS has built a catalog of loci and uses locus ID and position within the locus (stored in columns 1 and 6 in this file) to retrieve information on these SNPs (i.e. STACKS can't use chromosome and position to identify SNPs). Now that we have created a 'whitelist', we can rerun populations and generate F-statistics and input files for these 100 highly differentiated SNPs only.
```
cd ~/stacks
mkdir 
populations -t 16 -P ~/stacks/gstacks/ -M ~/scripts/popmap_canada.txt -O populations_canada_random_white --fstats --vcf --genepop --structure -W ~/stacks/populations_canada_random/high_fst.whitelist.tsv
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

#### With VCFtools
To calculate F<sub>ST</sub> in VCFtools you need to provide the list of samples for each population to include in the calculation. For the redeuced dataset, including only 80 individuals from Greenland and Canada, you will find the list of samples for each of the two lineages in the ```~/scripts``` folder
```
vcftools --vcf populations_2lin_random/populations.snps.vcf --weir-fst-pop ~/scripts/pop_canada.txt --weir-fst-pop ~/scripts/pop_greenland.txt  --out fst_2lin
```
If you have a reference genome and panel of SNPs dense enough (or whole genome resequencing data) you can calculate window estimates with 
```
vcftools --vcf populations_2lin_random/populations.snps.filter.vcf --weir-fst-pop ~/scripts/pop_canada.txt --weir-fst-pop ~/scripts/pop_greenland.txt  --out fst_2lin_win --fst-window-size 100000 --fst-window-step 100000
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
The calculation is very fast, and in addition to per-site estimates you get mean and weighted F<sub>ST</sub> estimates in the stdout, which you could save in a separate file by adding ```> fst_2lin.log``` at the end of the previous command.


-> fix faststructure
-> finish DAPC with different datasets


```
populations -V populations.snps.filter.vcf -M ~/scripts/popmap_2lin.txt -O ../populations_2lin_random_white --fstats --vcf --genepop --structure

mkdir ../populations_all_random_white/
populations -V populations.snps.filter.vcf -M ~/scripts/popmap_all.txt -O ../populations_all_random_white --fstats --vcf --genepop --structure

vcftools --vcf populations_2lin_random/populations.snps.filter.vcf --weir-fst-pop ~/scripts/pop_canada.txt --weir-fst-pop ~/scripts/pop_greenland.txt  --out fst_2lin
vcftools --vcf populations_2lin_random/populations.snps.filter.vcf --weir-fst-pop ~/scripts/pop_canada.txt --weir-fst-pop ~/scripts/pop_greenland.txt  --out fst_2lin_win --fst-window-size 10000 --fst-window-step 10000
vcftools --vcf populations_2lin_random/populations.snps.filter.vcf --site-pi --keep ~/scripts/pop_canada.txt --out pi_canada
vcftools --vcf populations_2lin_random/populations.snps.filter.vcf --site-pi --keep ~/scripts/pop_greenland.txt --out pi_greenland
```
### Structure with faststructure
STRUCTURE has been the most popular software for detecting and describing population structure since it was published 20 years ago. However, it was developed before we started genotyping thousands of markers at a time and it is therefore slow for large genomic datasets. fastSTRUCTURE (Raj eet al. 2014, Genetics) provides a much faster alternative for these datasets while it is based on a similar Bayesian framework. 

We'll use the ```.structure``` file created by ```populations```. The first 6 columns of the file will be ignored; these typically would include IDs, metadata, etc. so we need to have 6 dummy columns before the genotypes for each individual/locus. Note that this software only handles bi-allelic loci. The two alleles at each locus can be encoded as desired; however, missing data should be encoded as '-9' (also check https://rajanil.github.io/fastStructure/ for details). Below is a little script to prepare the file and run fastSTRUCTURE. 
```
cd populations_2lin_random_white *** populations_all_random_white *** populations_canada_random_white
tail -n +3 populations.snps.filter.p.structure > populations.filter.structure_nohead ### deletes first two rows
cut -f1-2 populations.structure_nohead > twocolumns
paste twocolumns twocolumns populations.structure_nohead | sed 's/\t0/\t-9/g' > populations.faststructure.str
python structure.py -K 2 --input=populations.faststructure --output=capelin.all --format=str 
```
Now you can modify the ```-K``` parameters and test different different numbers of populations.
With this script you can test what number of populations/cluster best describes your datset.
```
python chooseK.py --input=test/testoutput_simple
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
```

### DAPC
A Discrimant Analysis of Principal Components is a multivariate approach that merges a Pricipal Component Analysis (PCA) and a Discriminant Analysis (DA). A PCA aims to summarize the variation among individuals and it runs very fast, also on large datasets. However, it is not powerful at discriminating groups because it doesn't use any a priori information on grouping and intra-group variation can overwhelm inter-group variation. A DA, on the other hand, tries to summarize the variation among groups, while minimizing the variatin within groups. Thus, a DAPC takes the best of the two analyses to 

Download file ```populations_all_random_white/populations.snps.filter.p.structure``` onto your desktop. Change extension ```.structure```to ```.str``` and open R Studio.
```
### load adegenet package
library(adegenet)
### set working directory
setwd("~/Desktop")
### load dataset and convert it from structure to genind format
all_random<-import2genind("populations.snps.filter.p.str")
### This command will prompt questions about the structure of the files
```
![convert4dapc](https://github.com/clairemerot/physalia_adaptation_course/blob/master/images_tutorial/adegenet_loaddata.png)

Now we can run the DAPC
```
dapc_all<-dapc(all_random)
scatter(dapc_all)
compoplot(dapc_all)
```
You can change pretty much everything in your plot and the authors of adegenet have put together a great tutorial <https://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf>
