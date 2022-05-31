#Libraries
library(dplyr)
library(magrittr)
library(tibble)
library(ggplot2)
library(reshape2)

popmap <- read.table("info_samples.csv", h=T,sep=';')

#2. load geno data for the 12 canadian pop
geno.012_can <- read.table("pca/populations_can_random.012")[,-1] #load genotype matrix
geno.012_can.pos <- read.table("pca/populations_can_random.012.pos") %>% #load SNPs info
  mutate(., locus=paste(V1,V2,sep='_')) #create a new column for SNP info name (CHR + position)
geno.012_can.indv <- read.table("pca/populations_can_random.012.indv") #load individuals info

Set rownames and colnames to the geno matrix
dimnames(geno.012_can) <- list(geno.012_can.indv$V1, geno.012_can.pos$locus)
#check the geno matrix
geno.012_can[1:12,1:9]

geno.012_can[geno.012_can == -1] <- NA

geno.012_can.imp <- apply(geno.012_can,2,function(x){
  replace(x, is.na(x), as.numeric(names(which.max(table(x))))) })

geno.012_can.imp[1:12,1:9]

pca.can <- prcomp(geno.012_can.imp)
screeplot(pca.can)

#get stats info from the pca
sum.pca <- summary(pca.can)
#print stats info
sum.pca$importance[,1:5]

#prepare dataset to plot PCAs
pca.can.sub <- pca.can$x[,1:4] %>% #retain the first four PCs
  as.data.frame(.) %>% #transform to dataframe object
  tibble::rownames_to_column(., var='id') %>% #set rownames to a new column for samples ids
  dplyr::left_join(., popmap, by='id') #Here we use the left_join function
#from dplyr to wrap the population vector
#of our samples.

#plot by pop
ggplot(pca.can.sub) + aes(x=PC1, y=PC2, col=pop) +
  geom_hline(yintercept = 0, lty=2, col='grey50') + #add horiz line at y=0
  geom_vline(xintercept = 0, lty=2,col='grey50') +  #add vertical line at x=0
  geom_point() + 
  theme_bw() + #use  classic dark-on-light ggplot2 theme
  theme(panel.background = element_rect(fill='white'), #set some theme options
        panel.grid = element_blank())


#plot by sex
ggplot(pca.can.sub) + aes(x=PC1, y=PC2, col=sex) +
  geom_hline(yintercept = 0, lty=2, col='grey50') + #add horiz line at y=0
  geom_vline(xintercept = 0, lty=2,col='grey50') +  #add vertical line at x=0
  geom_point() + 
  theme_bw() + #use  classic dark-on-light ggplot2 theme
  theme(panel.background = element_rect(fill='white'), #set some theme options
        panel.grid = element_blank())

#look at lodings
loadings.melt <- reshape2::melt(abs(pca.can$rotation[,1:4])) %>% #get absolu teloadings values
  set_colnames(., c('SNPs','PC','loading')) %>% #set the colnames of the new dataframe
  mutate(., CHR=substr(SNPs,1,4)) #create a new column to inform about chromosome

#plot the data
ggplot(data=loadings.melt) +
  geom_bar(aes(x=SNPs, y=loading, fill=CHR), stat='identity') +
  facet_wrap(~PC) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
