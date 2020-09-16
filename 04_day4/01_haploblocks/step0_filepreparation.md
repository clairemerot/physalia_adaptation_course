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
