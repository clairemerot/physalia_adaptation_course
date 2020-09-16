## Step 1 Genotyping the individuals
This is the portion I cut which explains how I classified the individuals into three clusters


Thanks to our local PCA exploration, we know that there are non-recombining haploblocks which may be an inversion on chromosome 4.
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
