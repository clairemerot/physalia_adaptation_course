## Optional: Step4 Annotation of repeated regions and CNV
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

