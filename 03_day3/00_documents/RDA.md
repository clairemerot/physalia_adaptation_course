#### Prepare file (if you ahve not done so in the main tutorial)
```
geno<-read.table("02_data/geno_matrix.txt")
SNP_pos<-read.table("02_data/SNP_pos.txt", header=T)

#transpose data and give meaningful colnames
gen<-t(geno)
colnames(gen)<-paste(SNP_pos$chromosome, SNP_pos$position, sep="_")
gen [1:10,1:10]

#replace 9 by NA
gen[which(gen=="9")]<- NA
#evaluate % of missing
sum(is.na(gen))/(dim(gen)[1]*dim(gen)[2]) # about 1.5% of missing data
#impute missing with the most common geno
gen.imp <- apply(gen, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
sum(is.na(gen.imp)) # No NAs

#env data
info<-read.table("02_data/info_samples_canada.txt", header=T)
head(info)
```

####  environmental variation vs.geographic variation?
With the RDA we can also ask about the relative contribution of geography and environment since the rda can be controlled by a third matrix. It is also possible to include a previously inferred neutral structure, etc.

Very simply we can use latitude and longitude. Run the rda with this geogrpahic variables. As we have two RDA now, it may be worth looking at the triplot and explore significance axis by axis

```
geo.rda <- rda(gen.imp ~ info$lat + info$long, scale=T)
geo.rda
RsquareAdj(geo.rda)

jpeg("06_rda/rda1_triplot_geo.jpeg")
plot(geo.rda, scaling=3) 
points(geo.rda, display="sites", pch=20, cex=1.3, col=info$pop, scaling=3)
dev.off()
```
![rda1_triplot_geo](https://github.com/clairemerot/physalia_adaptation_course/blob/2022/03_day3/07_img_readme/rda1_triplot_geo.jpeg)

Only RDA1 describe significant association between genetic variation and latitude. Yet, geography explains only 0.2% of variance in our case!
Well, let's try the controlled the rda by geography anyway:

```
temp.geo.rda <- rda(X= gen.imp, Y= info$temp, Z=cbind(info$lat,info$long), scale=T)
temp.geo.rda
RsquareAdj(temp.rda)
RsquareAdj(temp.geo.rda)
```
So here, we can interpret that temperature explain 0.1% of variance, of which 0.08% is also explained by geography. I guess that in some systems, or with more populations, those numbers would be higher!!

#### Advanced options for geographic variables
Another option is to use dbmem, which generate vectors transforming (spatial) distances to rectangular data that is suitable     for constrained ordination or regression.
Since we have multiple variable (multiple pcnm vectors), we want to know which ones to keep. The function ordistep will perform model choice to select variable. It is a bit longer to run. 

```
library(cluster)
#Create a geographic matrix with dbmem
coorgeo<-info[,6:7]
dist_eucl<-daisy(coorgeo)
geo_pcnm<-pcnm(dist_eucl)
geo_mat<-geo_pcnm$vectors
head(geo_mat)

#perform rda
geo.rda<-rda(gen.imp~geo_mat[,1]+geo_mat[,2]+geo_mat[,3]+geo_mat[,4]+geo_mat[,5]+geo_mat[,6]+geo_mat[,7]+geo_mat[,8])
RsquareAdj(geo.rda)
#select model
anova(geo.rda, step=1000, by= "margin")
ordistep(geo.rda)
```
We now explain more than 1.15% of variance with geographic distance... 

#### Repeat the analysis on phenotype?
If you are interested in finding SNPs/multiloci haplotypes associated with phenotype, this can also be done with rda.
Here, you can pratice with sex, and check whether the SNPs are indeed located on Chr5!
```
sex.rda <- rda(gen.imp ~ info$sex, scale=T)

RsquareAdj(sex.rda)

sex.signif.full <- anova.cca(sex.rda, parallel=getOption("mc.cores")) # default is permutation=999
sex.signif.full

load.sex.rda <- scores(sex.rda, choices=c(1), display="species") 
load.sex.rda.pos<-cbind(load.sex.rda, SNP_pos)

#plot
jpeg("06_rda/loading_sex_manhattanplot.jpeg")
ggplot(load.sex.rda.pos, aes(x=position, y=RDA1, colour=chromosome))+ 
  geom_point()+
  theme_classic()+
  facet_grid(cols = vars(chromosome), scales = "free_x", space="free_x")
dev.off()
```

#### expand to multivariate?
Another powerful aspect of the RDA is to be able to analyse together several environmental variables. 

Yet, they should not be correlated. so if you include several environmental axis, it is best to first do quick correlation steps and select the most interesting variables or PC axis that summarize them. 

The problem of including several variables (or worse sort of composite variables based on PCA is that it can become tricky to disentangle which outliers are associated with which environmental variable). 

Today, we only included temperature but if you wish, you can use the gps points and databases (bioclim, marspec...) to retrieve more environmental data. If you are interested into that we can come back on it later or on the last day. We palso resent here a sort of dummy example with latitude and longitude as variables if you have time to run the multivariate.

If you include several environmental variables (for instance rda(gen.imp~info$temperature +info$longitude), you can test of teh association between genetic variation and each variable with the model by axis, and check that the correlation between environmental variable with the VIF. It might take a bit longer to run. The ordistep function shown above will be useful to select meaningful variables.

To extract outliers associated with each environmental variable, you can have a look at Forester's tutorial.

```
temp_gps.rda <- rda(gen.imp ~ info$temperature+ info$lat + info$long, scale=T)
temp_gps.rda

RsquareAdj(temp_gps.rda)

signif.axis <- anova.cca(temp_gps.rda, by="axis", parallel=getOption("mc.cores"))
signif.axis
vif.cca(temp_gps.rda)
```
