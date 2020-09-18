Today, we will explore the annotation of the genome, to see if some of the outliers we found might belong to candidate genes, the distribution of our SNPs in exons, introns, regulatory regions,etc. We will also look at the annotation of repeated and transposable elements and compare that to our CNV detection.

## Getting started 

First of all, like previous days please copy into your directory the folder 05_day5, open it and we will work from there all day.

```
~/Share/physalia_adaptation_course/05_day5 .
cd 05_day5
```

Inhere, you should have everything we need for today. You can explore the different folders and files
``` ls 02_data``` in which you will find the vcf, ``` ls 03_outliers``` in which I put the txt files that we exported on day 3 when analysing the association with temperature, and on day4 when analysing divergence between haploblocks or between sexes
I also put here the list of all SNPs present in the vcf. you can look at each file with the command ``` head 03_outliers/SNP_pos.txt``` for instance

The annotated transcriptome (which is generally in .gff format) as well as the transposition onto the genome in an annotation table are located at the following path ~/Share/ressources/
You don't need to copy them as we have prepared simplified files (in R simply selecting the relevant column) when needed. You may want to have a loof to get a sense of what it looks like.
``` less ~/Share/ressources/genome_mallotus_dummy.gff3``` 
press "q" to exit the less visualisation
``` less ~/Share/ressources/genome_mallotus_dummy_annotation_table.tsv``` 

## Step 1 SNPeff : annotating our snps 

SNPeff is a program that uses the gff file and the position of each SNP to annotate the vcf.
If you work on a model species which already has a database, you are lucky! If not, you need to build a database. 

As it is a bit long, and takes space on the server, I have done it for you. It was not very straightforward, so I thought I will keep track of how I did and put it for you in this file https://github.com/clairemerot/physalia_adaptation_course/blob/master/05_day5/SNPeff_createDB (but do not try to run it on the server please)

If you want to, you can look at the database by doing:

```
java -jar ~/Share/ressources/snpEff/snpEff.jar dump genome_mallotus_dummy | less
```
To exit "less" simply press "q"

#### Annotate the vcf 
Now we can annotate our vcf. As you are getting used to now, we use a raw vcf file in the folder 02_data and will write the output into the folder 04_snpEff in which we will have all subsequent files related to the snpEff analyses
```
java -Xmx4g -jar ~/Share/ressources/snpEff/snpEff.jar genome_mallotus_dummy 02_data/canada.vcf > 04_snpEff/canada_annotated.vcf
```
Let's look at the new vcf
```
less -S 04_snpEff/canada_annotated.vcf
```
As you see it kepts the vcf format with its specific header but this is not easy to use as it is now if we want to import SNP information in R.
We will use a few bash command and awk to spli the information by the symbol "|" and rather make different column separated by tab

```
cat 04_snpEff/canada_annotated.vcf | awk -F"|" '$1=$1' OFS="\t" | cut -f 1-9 > 04_snpEff/SNP_annotated_formatted.txt
#look at the first 25 lines
head -n 25 04_snpEff/SNP_annotated_formatted.txt
#look at the last lines
tail 04_snpEff/SNP_annotated_formatted.txt
```

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
we can join both matrix of repartition to have a look:
```
joined_repartition1<-as.data.frame(cbind(as.matrix(table(snpEff_db$category)),as.matrix(table(outlier_annotated$category))))
colnames(joined_repartition1)<-c("n_all_snps","n_oulier")
joined_repartition1$category<-as.character(row.names(joined_repartition1))
joined_repartition1
```
And we can use Fisher exact test to ask whether one category of variant is enriched in our outlires
```
total_snp<-dim(snpEff_db)[1]
total_outlier<-dim(outlier_annotated)[1]

# for row i (category i)
i=1
joined_repartition1$category[i]

#build the contigency table
cont_tab<-rbind(c(joined_repartition1$n_oulier[i],total_outlier),
                c(joined_repartition1$n_all_snps[i],total_snp))
cont_tab

#run Fisher test
cont_tab[ which(is.na(cont_tab))]<-0
fish_res<-fisher.test(cont_tab,alternative = "greater")
fish_res
```
You can change Fisher alternative hypothesis to "two.sided", "greater" or "less"

You can loop over the categories with a loop and store the results
```
for (i in 1 : length(joined_repartition1$category)
{
cat(joined_repartition1$category[i])

#build the contigency table
cont_tab<-rbind(c(joined_repartition1$n_oulier[i],total_outlier),
                c(joined_repartition1$n_all_snps[i],total_snp))
cat(cont_tab)

#run Fisher test
cont_tab[ which(is.na(cont_tab))]<-0
fish_res<-fisher.test(cont_tab,alternative = "greater")
cat(fish_res)
}
```



## Step 2 Bedtools : find the intersection between SNPs and genes
Bedtools is a program that is great to find the intersection between two files. It usually works on a specific bedformat which have at least three columns (Chromosome, FromPosition, ToPosition) and 12 columns to the maximum.
It won't like having header. We will try to keeep a 4th column with SNP id.

I have prepared for your the annotation file of the genome in a bed-readable format, that you may want to see by doing
```
less -S 05_bed/genome_mallotus_dummy_annotation_simplified.bed
```
(q to exit from less)
q
If you remember our outliers files they were not formatted as such... Because we do not cover all the genome, we will look for genes in a window of X kb around the SNP position. The size of this window should ideally be adjust depending on LD decay in your organism (which can be assessed by plotting LD against distance with the output of plink that we found yesterday... For today we will choos 10 kb but if you are curious you cna explore different size.
To prepare the files, we will use R. 

### In R (either on Rstudio in your computer or in a R terminal on the server)

```
#load file
outlier_temp_rda<-read.table("03_outliers/outlier_temp_rda.txt", header=T)
#have a quick look
head(outlier_temp_rda)
#what's its dimension?
dim(outlier_temp_rda)

#which size around the SNP
window<-10000
#add a vector with start position
outlier_temp_rda$start<-outlier_temp_rda$position-(window/2)
#oups it can't be ngative! replace negative by 0
outlier_temp_rda$start[outlier_temp_rda$start<0]<-0 
#add a vector with stop position
outlier_temp_rda$stop<-outlier_temp_rda$position+(window/2)
#have a look
head(outlier_temp_rda)

#which columns shoud we keep?
outlier_temp_rda_bed<-outlier_temp_rda[,c(2,5,6,1)]
#save your file
write.table(outlier_temp_rda_bed, "05_bed/outlier_temp_rda.bed", row.names=F, sep="\t", quote=F,col.names=F)
```
Now that we are good for the rda outliers, we also need to do so for all the SNPs if we wan tto do later gene enrichment analysis. We want to know the pool of genes in which outliers may have been drawn.

```
all_snps<-read.table("03_outliers/SNP_pos.txt", header=T)
all_snps$start<-all_snps$position-(window/2)
all_snps$start[all_snps$start<0]<-0 
all_snps$stop<-all_snps$position+(window/2)
head(all_snps)
all_snps_bed<-all_snps[,c(2,4,5,1)]
write.table(all_snps_bed, "05_bed/all_snps.bed", row.names=F, sep="\t", quote=F,col.names=F)
```

And you can do the same for the outliers from baypass, or the outliers of divergence on chr4 or chr5... (if you are late skip)

### on the server
If you did the bedfiles in Rstudio, Please copy back your formatted outliers and snp bedfiles into the 05_day5/05_bed folder on the server. 
We will now run bedtools. The command intersect will look for overlap between the file givne with -a, and the file given with -b, the argument -wb will print the information coming from the annotation file (gene names, gene ontology, uniprot ID, etc)
the command ">" redirect the output in the file of your choice


With the command wc -l we will count the lines to see how many transcripts intersect with our outliers

```
bedtools intersect -a 05_bed/outlier_temp_rda.bed -b 05_bed/genome_mallotus_dummy_annotation_simplified.bed -wb > 05_bed/outlier_temp_rda.intersect
cat 05_bed/outlier_temp_rda.intersect | wc -l
head 05_bed/outlier_temp_rda.intersect
less -S 05_bed/outlier_temp_rda.intersect
```
For gene enrichment analyses we will need a simpler format, so we will cut the 8th column (transcript name) with the command cut and keep it for goseq
```
cat 05_bed/outlier_temp_rda.intersect | cut -f 8 > 06_go/outlier_temp_rda.transcript
head 06_go/outlier_temp_rda.transcript
```

We repeat the analyses for all SNPs

```
bedtools intersect -a 05_bed/all_snps.bed -b 05_bed/genome_mallotus_dummy_annotation_simplified.bed -wb > 05_bed/all_snps.intersect
cat 05_bed/all_snps.intersect | wc -l

cat 05_bed/all_snps.intersect | cut -f 8 > 06_go/all_snps.transcript
head 06_go/all_snps.transcript
```

You may now repeat the commands for other outliers files if you wish. (if you are late skip)

### Have a look at your results
So copying back your .intersect files (and .transcript) on your computer (inside 05_day5/05_bed of course!), you can easily open them with a text editor (notepad ++ or equivalent). You cna see the name of the genes, gene ontology. What do you think?
For information about gene ontology codes, you may want to check this website http://geneontology.org/docs/ontology-documentation/

For information about proteins, you can look at https://www.uniprot.org/

## Step3 Goseq Gene ontology enrichment analysis
#### WARNING #### Gene enrichment analysis are mostly designed for RNAseq analyses in which all genes are analysed and of much better use when you have whole-genome data. Here on RAD-seq data, in which we subsample the genome, and with so few outliers, we should not make this kind of analysis, or at least not overinterpret it. We decided to include it in the tutorial for teaching purpose since most of you will be working with whole genome data but please keep in mind this warning!

Note: This script was built with the help or Dr. E. Berdan (U. Stockholm)

This part will only be in Rstudio on your computer.
Please install the library goseq ( https://bioconductor.org/packages/release/bioc/html/goseq.html)
Disclaimer: I admit that this library is a pain because you will see that it requires many little twists like putting the transcripts in row.names, etc that makes it difficult to prepare the data... Hopefully there will be simpler versions or maybe one of you knows an easier tool that he/she may recommend...


```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("goseq")
library(goseq)
```


####  Prepare our correspondance transcript/GO
Again, as we are on a non-model organism (and a dummy genome), there is no database of our transcriptome and GO terms in the library, so we will build it.

We have prepared for you a simplified annotation of the transcripts, and tried to remove duplicates. In fact, when a genome gets annotated, or when contigs from a transcriptome are aligned to a genome, one may have mutliple matches, which in turn may inflate our calculations. (if one SNP hits 10 times the same gene but with a different transcript name..)

```
transcript_info = read.table("06_go/transcript_go_simplified.txt", header=T, stringsAsFactors=F, sep="\t")
head(transcript_info)
row.names(transcript_info)<-transcript_info$TranscriptName
```
As you see go terms are all side by side, which will not be super helpful to make a list transcripts/GO
Please install the library splitstackshape which is able to split a whole column, and the library data.table which is super helpful to transform wide datatbale into long datatable
```
library(splitstackshape)
#split the GO terms
go_split = cSplit(transcript_info, "GeneGo", sep = ";")
go_split$contig<- transcript_info$TranscriptName
head(go_split)
```

```
library(data.table)
#linearize the matrix
terms = colnames(select(go_split, contains("GeneGo")))
go_long = melt(go_split,measure.vars=terms, id.vars = "contig", na.rm=TRUE)
head(go_long)
go_ready = as.data.frame(go_long[,c(1, 3)])
head(go_ready)
```

####  Enrichment test
First we need to know which genes were in our pool of analysis, in other words, which genes are nearby a SNP covered by a RAD loci. This will be a subset of all the genes presents in the transcriptome/genome since we did a reduced-representation sequencing. Note that for whole genome you may be in the same situation if because of coverage, quality or whatever reason we have genes not covered by SNPs (after filtration)

```
#upload transcript intersecting with snps
all_transcripts<-read.table("05_bed/all_snps.transcript", header=F)
colnames(all_transcripts)[1]<-"TranscriptName"
dim(all_transcripts) #how many?
head(all_transcripts)
```

Now we want to add the information about the size of the gene, as a way to correct for bias that one is more likely to fall inside a long gene than a short one .
We will use the left_join command from dplyr library which is super handy. It will match by transcript name and add the rest of the table, matching rows (magi function to remmebmber!)

```
library(dplyr)
#add size info
all_transcripts<-left_join(all_transcripts,transcript_info[,c(1,3)])
dim(all_transcripts)
head(all_transcripts)
```
Now a problem is that some transcripts are listed twice, it could be that the same gene matches at several places in the genome, that we have two snps in the same gene, etc. So we will use a finction to remove duplicated rows. (magic function of dplyr again!)
```
#make unique
all_transcripts_unique<- all_transcripts %>% distinct(TranscriptName,.keep_all = TRUE)
dim(all_transcripts_unique) # see the matrix has reduced...
#and we need to name the rows b transcript names for goseq...
row.names(all_transcripts_unique)<-all_transcripts_unique$TranscriptName
head(all_transcripts_unique)
```
Now let's open one of our outlier list and format it. We will start with the outliers from the associations with temperature as found by the rda

```
#transcripts in outliers
outliers_transcripts_temp_rda<-read.table("05_bed/outlier_temp_rda.transcript", header=F)
colnames(outliers_transcripts_temp_rda)<-"TranscriptName"
head(outliers_transcripts_temp_rda)
dim(outliers_transcripts_temp_rda)
```
We will now add a column to the matrix with all the genes indicating whether this gene is found or not inside the outliers list. this will be a 0/1 vector
```
all_transcripts_unique$outliers_temp_rda<- as.numeric(all_transcripts_unique$TranscriptName %in% outliers_transcripts_temp_rda$TranscriptName)
head(all_transcripts_unique)
```

Now we run goseq function (nullp) to prepare the data and integrate length bias. It requires vectors so here we go transforming data:

```
measured_genes = as.vector(all_transcripts_unique$TranscriptName)
outliers_genes = as.vector(all_transcripts_unique$outliers_temp_rda)
length = as.vector(all_transcripts_unique$length)
pwf_outliers = nullp(outliers_genes, bias.data=length)
row.names(pwf_outliers)<-row.names(all_transcripts_unique)
head(pwf_outliers)
```
Yeah! we have formatted all the files so now, let's test enrichment using our database (transcript/GO) named "go_ready" and the prepared list of genes with 0/1 info for our outliers from temperature association with rda 
```
enrich_outliers = goseq(pwf_outliers,gene2cat=go_ready, use_genes_without_cat=TRUE)
head(enrich_outliers)
```
We will now correct for multiple testing with Benjamini & Hochberg correction
```
enrich_outliers$over_represented_padjust<-p.adjust(enrich_outliers$over_represented_pvalue,method="BH")
head(enrich_outliers)
write.table(enrich_outliers, "06_go/GO_enrich_temp_RDA.txt", sep="\t")
```
all go terms are presented in this matrix.We ahev exported it so you can look at it in an editor.

But we may want to see only the the significant enrichments:
```
enrich_outliers[which(enrich_outliers$over_represented_padjust<0.10),]
```
What do you think? Sometimes it may be significan tbut when one has low numbers (1 gene out of 4), is this really interpretable? 
Anyhow, you know how to do it!

If you wish, you can repeat the analysis on the outliers from baypass, or joined BP/rda, or on the outliers of divergence found on chr4 and chr5. 

## Step4 Annotation of repeated regions and CNV
I have run a programm called RepeatMasker http://www.repeatmasker.org/ which uses databases of transposable elements and detection of repeated pattern to annotate the genome for those interspersed repeats and low-complexity DNA sequences. Because we are on afish, I use the database Danio rerio (Zebrafish). The command is very easy
[do not run] ```RepeatMasker genome_mallotus_dummy.fasta -pa $N_CPU -species Danio```

Another possibility is to build the database of TE and repeats for your species. A good programm for that is RepeatModeler2 (Flynn et al PNAS 2020 https://doi.org/10.1073/pnas.1921046117)

So RepeatMasker outputs both a masked version of the genome and a file .out in which there is a list of repeated elements and TE with their position. I have arranged a little that file to make it into a bed format. you can see it in the folder 07_TE
``` head 07_TE/genome_mallotus_dummy_repeat.bed ```

Now we will look into the loci that we identified as being duplicated (CNVs). Yann puts back for you the CNV matrix in the 07_TE folder so don't bother looking at yesterday's files :-)

We will use R to convert the full matrix and the matrix of outliers into bed format. Please open R either in the terminal or  on your own computer (but you will need to move files between server and and local).
We will choose a window of 500 bp around the CNV locus to see if this locus may belong to a repeated element.

```
#load data
 CNV<-read.table("07_TE/caplelin_NWA_list_of_CNVs_loci.txt", header=T)
 head(CNV)
 
 window=500
 #calculate start adn stop of the window around the CNV loci
 CNV$start<-CNV$POS-(window/2)
 CNV$start[CNV$start<0]<-0
 CNV$stop<-CNV$POS+(window/2)
 
 #subset the matrix of outliers
 CNV_outlier<-CNV[CNV$outlier=="TRUE",]
 head(CNV_outlier)

#format for bedtools
CNV.bed<-CNV[,c(2,5,6,1)]
CNV_outlier.bed<-CNV_outlier[,c(2,5,6,1)]
write.table(CNV.bed, "07_TE/CNV.bed", row.names=F, col.names=F, sep="\t", quote=F)
write.table(CNV_outlier.bed, "07_TE/CNV_outlier.bed", row.names=F, col.names=F, sep="\t", quote=F)
```

Now we can run bedtools to output the intersection of the CNV positions with the annotated repeats
```
bedtools intersect -a 07_TE/CNV.bed  -b 07_TE/genome_mallotus_dummy_repeat.bed -wb > 07_TE/CNV_repeat.intersect
cat 07_TE/CNV_repeat.intersect | wc -l
head 07_TE/CNV_repeat.intersect

bedtools intersect -a 07_TE/CNV_outlier.bed  -b 07_TE/genome_mallotus_dummy_repeat.bed -wb > 07_TE/CNV_outlier_repeat.intersect
cat 07_TE/CNV_outlier_repeat.intersect | wc -l
cat 07_TE/CNV_outlier_repeat.intersect
```

And you can now have a look at your results with a text editor and investigate the composition of your detected loci in different class of TE with R. 

As you see, maybe our window of 1000bp was a bit wide, sometimes the CNV intersect with several simple_repeat regions... Maybe this is also a sign that some of them occur in this kind of area? 

We could also have looked at the intersect between CNV and genes by running bedtools with the annotated genome as we did earlier.

