# Date: Wed Aug 26 16:54:06 2020
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
  library(vegan)
# --------------

# Set working dir -------------------------
#setwd(PATH)

##########################################################################################
######################        Define global functions       ##############################
##########################################################################################

# load specific useful functions from source
'%ni%' <- Negate('%in%') # reverse of %in%

#---------- Add new project functions -------------------------------

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

##########################################################################################
###########################        Main script       #####################################
##########################################################################################


# |---------|
# | Step 1  | ================> load data
# |---------|

#1.load population map
popmap <- read.table("01-info_files/info_samples.csv", h=T, sep=';') %>%
  select(., -file) %>% #remove column file -> unecessary
  filter(., region == 'NWA')

#2.import CNVs matrix format
CNVs.mat <- read.table("02-data/capelin_CNVs_gdepth_normalized.txt") %>%
  t(.) # transpose matrix indv x loci
CNVs.mat[1:4,1:4]

#remove loci from the sex chromosome
CNVs_no_sex_chr.mat <- CNVs.mat[,!grepl('Chr5_', colnames(CNVs.mat))]

#Impute missing values with the average of norm. read count overall
CNVs_no_sex_chr.imp.mat <- apply(CNVs_no_sex_chr.mat,2,function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
CNVs_no_sex_chr.imp.mat[1:4,1:4]

# |---------|
# | Step 2  | ================> Make RDA
# |---------|

#perform RDA
rda <- vegan::rda(CNVs_no_sex_chr.imp.mat ~ popmap$temperature, scale=T)

#check RDA R2 adjusted
vegan::RsquareAdj(rda)
#check RDA significance
anova(rda)
vegan::anova.cca(rda, by="axis",permutations = 1000, parallel = 6)

# |---------|
# | Step 3  | ================> Find outliers from RDA
# |---------|

load.rda <- scores(rda, choices=c(1), display="species")  # Species scores for the first constrained axes

cand <- outliers(load.rda[,1],2.25) #required function to find outliers
ncand <- length(cand)

cand <- cbind.data.frame(rep(1,times=length(cand)), names(cand), unname(cand))

colnames(cand) <-  c("axis","locus","loading")

cand$snp <- as.character(cand$snp)

pred <- dplyr::select(popmap, temperature)

foo <- matrix(nrow=(ncand), ncol=1)  # n columns for x predictors
colnames(foo) <- c("correlation")

#get corrleation level between CNVs distribution and environment
for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <-CNVs_no_sex_chr.imp.mat[,nam]
  foo[i,] <- apply(pred,2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,foo)
head(cand)

#define list object of candidate CNVs loci
list_of_CNV_candidates <- cand$locus
###### end RDA analysis


list_of_CNVs_loci <- data.frame(CNV_locus=colnames(CNVs_no_sex_chr.imp.mat)) %>%
  mutate(., CHR=gsub('_.*', '',CNV_locus),
            POS=gsub('.*_','',CNV_locus),
            outlier = ifelse(CNV_locus %in% list_of_CNV_candidates,TRUE,FALSE))

head(list_of_CNVs_loci)

# |---------|
# | Step 4  | ================> Check adaptive structure lead by outliers loci
# |---------|

#filter the initial CNV matrix
CNVs_no_sex_chr.outliers.mat.imp <- CNVs_no_sex_chr.imp.mat[,colnames(CNVs_no_sex_chr.imp.mat) %in% list_of_CNV_candidates]

# perform PCA with candidate CNVs associated with temperature
pca.outliers <- prcomp(CNVs_no_sex_chr.outliers.mat.imp, scale=T)

#check the variances against the number of the principal component. 
#png("screeplot_CNVs_outliers.png")
screeplot(pca.outliers)
#dev.off()

# Prepare dataframe for PCA
pca.outliers.sub <- pca.outliers$x[,1:3] %>%  #select only the three first PCs
  as.data.frame(.) %>%  #define matrix as dataframe
  rownames_to_column(., var='id') %>% #mutate rownames to unique column of samples ids
  mutate(., pop = popmap$pop, # add pop info
            temperature=popmap$temperature) %>% #add temperature info
  arrange(., desc(temperature)) %>%  #sort temperatures decreasing order
  mutate(., color=rep(colorRampPalette(c('cyan',"darkblue"))(12),rep(20,12))) #add color scale 

#Get sites centroids positions for plot
sites_centroids <- stats::aggregate(.~pop,data = pca.outliers.sub[,c(2,3,5)], mean) %>% 
  left_join(., select(pca.outliers.sub, pop, color), 'pop') %>% #add color info
  distinct(.) #get only sites info

#check centroids dataframe
sites_centroids 

#plot the results
ggplot() + 
  geom_hline(yintercept = 0, lty=2, col='grey50') +
  geom_vline(xintercept = 0, lty=2,col='grey50') +
  geom_point(data=pca.outliers.sub, aes(x=PC1,y=PC2, col=pop),cex=1)  + #add samples points
  scale_color_manual(values=as.character(as.factor(sites_centroids$color))) + #set colors
  geom_point(data=sites_centroids,aes(x=PC1,y=PC2,fill=pop),cex=3, pch=23) + # add sites centroids
  scale_fill_manual(values=sites_centroids$color) + #set colors
  geom_text(data=sites_centroids ,aes(x=PC1,y=PC2), label=sites_centroids$pop, #add site names
            nudge_x = 0.0, nudge_y = -0.20) +
  theme_bw()+ #set black&white ggplot theme 
  theme(panel.background = element_rect(fill = 'white', color='black', size = 1.2),
        panel.grid = element_blank(),
        axis.title = element_text(size=12, face='bold'),
        axis.text = element_text(size=12, face='bold'))

ggsave("PCA_biplot_oultiers_CNVs.png", width = 8,height = 6)
        
