# Day2: Population structure and differentiation 
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
