# Date: Fri Aug 28 2020
# --------------
# Author: Yann Dorant
# Date:
# Modification:
# --------------
# Libraries
library(dplyr)
library(magrittr)
library(tibble)
library(ggplot2)
library(reshape2)
# --------------

# Set working dir -------------------------
#setwd(PATH)

##########################################################################################
###########################        Main script       #####################################
##########################################################################################


# |---------|
# | Step 1  | ================> load data
# |---------|

#1. load population map
popmap <- read.table("01-info_files/info_samples.csv", h=T,sep=';')

#2. geno data 2lin.snps.vcf
geno.012_2lin.012 <- read.table('02-data/population.2lin.rand.snp.012')[,-1]
geno.012_2lin.012.pos <- read.table('02-data/population.2lin.rand.snp.012.pos') %>% 
  mutate(., locus=paste(V1,V2,sep='_'))
geno.012_2lin.012.indv <- read.table('02-data/population.2lin.rand.snp.012.indv')

#set rownames and colnames of the geno matrix
dimnames(geno.012_2lin.012) <- list(geno.012_2lin.012.indv$V1, geno.012_2lin.012.pos$locus)

#check the data matrix
geno.012_2lin.012[1:6,1:6]

# |---------|
# | Step 2  | ================> Impute missing data
# |---------|

# Here we will fill the NAs values by the most common genotype
# across all samples for a given SNP

# Missing data in 012.geno files from vcftools are coded by -1.
# We will change it for NAs
geno.012_2lin.012[geno.012_2lin.012 == -1] <- NA

#Then, impute missings data
geno.012_2lin.012.imp <- apply(geno.012_2lin.012,2,function(x){
  replace(x, is.na(x), as.numeric(names(which.max(table(x))))) })

#Remove ultra rare alleles i.e. SNPs alleles only represented by one individual
MAS <- function(geno, MAS.thresh =2){
  #get_MAS_locus <- apply(geno,2,function(x) length(which(x > 0)))
  get_MAS_locus <- apply(geno,2,function(x) sum(table(x)[-1]))
  blacklist <- get_MAS_locus[get_MAS_locus < MAS.thresh]
  message(length(blacklist)," SNPs removed")
  return(geno[,-which(colnames(geno) %in% names(blacklist))])
}

#apply MAS filter (remove singletons or monomorphic)
geno.012_2lin.012.imp.MAS2 <- MAS(geno.012_2lin.012.imp, 2)

# |---------|
# | Step 3  | ================> Perform PCA
# |---------|

#1 PCA line code
pca.2lin <- prcomp(geno.012_2lin.012.imp.MAS2, scale=T)

# #check the variances against the number of the principal component.
screeplot(pca.2lin)

summary(pca.2lin)$importance[,1:4]
#                           PC1      PC2      PC3      PC4
#Standard deviation     17.66092 7.271187 7.180096 7.033057
#Proportion of Variance  0.09911 0.016800 0.016380 0.015720
#Cumulative Proportion   0.09911 0.115910 0.132290 0.148010


#prepare dataset to plot PCAs
pca.2lin.sub <- pca.2lin$x[,1:4] %>% #retain the first four PCs
  as.data.frame(.) %>% #transform to dataframe object
  tibble::rownames_to_column(., var='id') %>% #set rownames to column samples ids
  dplyr::left_join(., popmap, by='id') #Here we use the left_join function
                                       #from dplyr to wrap the population vector
                                       #of our samples.

#Plot the PCA biplot
ggplot(pca.2lin.sub) + aes(x=PC1, y=PC2, col=pop) +
  geom_hline(yintercept = 0, lty=2, col='grey50') + #add horiz line at y=0
  geom_vline(xintercept = 0, lty=2,col='grey50') +  #add vertical line at x=0
  geom_point() + #add the samples
  scale_color_manual(values=c('#cab2d6','purple','royalblue', 'cyan')) + #set a new color scale
  theme_bw() + #use  classic dark-on-light ggplot2 theme
  theme(panel.background = element_rect(fill='white'), #set some theme options
        panel.grid = element_blank())

#Print it !
ggsave("PCA_biplpot_2lin_capelin.png", width = 6, height = 5)

# |---------|
# | Step 4  | ================> Expolore PCA loadings
# |---------|

#prepare dataset
loadings.melt <- reshape2::melt(abs(pca.2lin$rotation[,1:4])) %>% #get absolu teloadings values
  set_colnames(., c('SNPs','PC','loading')) %>% #set the colnames of the new dataframe
  mutate(., CHR=substr(SNPs,1,4)) #create a new column to inform about chromosome

#plot the data
ggplot(data=loadings.melt) +
  geom_bar(aes(x=SNPs, y=loading, fill=CHR), stat='identity') +
  facet_wrap(~PC) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

#Print it !
ggsave("PCA_loadings.png", height = 4, width = 15)
