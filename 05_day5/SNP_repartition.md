This tutorial makes use of SNPeff annotation

### On your computer, in R studio
Please copy *the whole folder 05_day5* to your local computer, and set in Rstudio your working directory as 05_day5.

We will quickly look into our snps annotation. 

First, let's open the whole list of SNPs and look how many introns, UTR, intergenic regions, etc we had in our dataset:

```
#load library
library(dplyr)

#open file
snpEff_db<-read.table("04_snpEff/SNP_annotated_formatted.txt", header=F) 
head(snpEff_db)

#simplify the data
snpEff_db<-snpEff_db[,c(1,2,3,9)] 

#put informative col names
colnames(snpEff_db)<-c("chr","position","id_snp","category")
head(snpEff_db)

#count how many of each variant category
#with Rbase
all_snps_repartition<-as.matrix(table(snpEff_db$category))
all_snps_repartition
```
We can now look at one of our outliers' list and use the magic function from dplyr to match snp_id and annotate our outliers
```
#load data
outlier_temp_rda<-read.table("03_outliers/outlier_temp_rda.txt", header=T)
head(outlier_temp_rda)

#use inner-join function to match the two datatable and keep only rows in common

outlier_annotated<-inner_join(outlier_temp_rda,snpEff_db)
head(outlier_annotated)

#count how many of each variant category
#with Rbase
outlier_repartition<-as.matrix(table(outlier_annotated$category))
outlier_repartition
```
