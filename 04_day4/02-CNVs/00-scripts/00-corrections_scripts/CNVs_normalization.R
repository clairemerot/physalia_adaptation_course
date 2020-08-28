# Date: Sept 01 2020
# --------------
# Author: Yann Dorant
# Modification:
# --------------
# Libraries
library(dplyr)
library(magrittr)
library(tibble)
library(edgeR)
# --------------

# Set working dir -------------------------
setwd(PATH)

##########################################################################################
######################        Define global functions       ##############################
##########################################################################################

# load specific useful function
'%ni%' <- Negate('%in%') # reverse of %in%

#---------- Add new project functions -------------------------------

##########################################################################################
###########################        Main script       #####################################
##########################################################################################

# |---------|
# | Step 1  | ================> load data files
# |---------|

#1.load population map
popmap <- read.table("01-info_files/info_samples.csv", h=T, sep=';') %>%
              select(., -file) #remove column file -> unecessary

CNVs_radTag <- read.table("CNVs_RADtag_IDs.txt")

#2. load gdepth file of duplicated SNPs (output from vcftools --geno-depth)
gdepth_raw_data <- read.table("02-data/capelin_NWA_4_60_0_6.gdepth", h=T) %>% 
                        mutate(., RadTag = paste(CHROM, CNVs_radTag$V1, sep="_")) %>%
                        select(., CHROM, POS, RadTag, everything())

#check the data
gdepth_raw_data[1:10,1:10]

#correct missing data -1 to 0 
gdepth_raw_data[gdepth_raw_data == -1] = 0

#check if missing value have been change to 0
gdepth_raw_data[1:10,1:10]

# |---------|
# | Step 2  | ================> prepare gdepth data for normalization
# |---------|

gdepth_transformed <- dplyr::select(gdepth_raw_data, -CHROM, -POS) %>% #remove POS col
  dplyr::distinct(., RadTag, .keep_all=TRUE) %>% #remove duplicated loci
  column_to_rownames(., var = 'RadTag')

# check data
dim(gdepth_transformed)
gdepth_transformed[1:10,1:10]

# |---------|
# | Step 3  | ================> perform normalization with EdgeR package
# |---------|
####################
DGE_list <- DGEList(counts=as.matrix(gdepth_transformed)) #create list to store info

gdepth_norm_Fact <-calcNormFactors(DGE_list)

gdepth_normalized <- cpm(gdepth_norm_Fact, normalized.lib.sizes = TRUE, log = F)

gdepth_normalized[1:10,1:10]

#change 0 values to NA
gdepth_normalized[gdepth_normalized == 0] <-NA

gdepth_normalized[1:10,1:10]

#write normalized matrix of CNV read depth
write.table(gdepth_normalized, "02-data/capelin_CNVs_gdepth_normalized.txt",
            col.names = T, row.names = T, quote = F, sep="\t")


#make tidy format of gdepth normalized read depths
gdepth.tidy <- reshape2::melt(gdepth_normalized) %>% set_colnames(., c("locus",'id','gdepth.norm')) %>%
  left_join(., popmap, by='id') %>% select(., id, pop, locus, gdepth.norm)
head(gdepth.tidy)

#write normalized dataframe of CNV read depth
write.table(gdepth.tidy, "02-data/capelin_CNVs_gdepth_normalized.tidy.txt",
            col.names = T, row.names = F, quote = F, sep="\t")
