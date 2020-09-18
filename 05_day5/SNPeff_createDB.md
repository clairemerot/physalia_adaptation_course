#### Prepare the database for SNPeff
First I have downloaded snpEff, unzip it and open it.
Then I have added a line into the config file, for this I usually use the editor nano which is easy-to-use.
You may see 
```
#new genome
genome_mallotus_dummy.genome : capelin
```
if you do 
```
less ~/Share/ressources/snpEff/snpEff.config
```
I putted the transcriptome called "genes.gff" into a data folder. and the genome as a fasta (.fa and .genome) in a genomes folder.
```
cp ~MYPATH/genome_mallotus_dummy.gff data/genome_mallotus_dummy/genes.gff
cp ~MYPATH/genome_mallotus_dummy.fasta data/genomes/genome_mallotus_dummy.fa
cp ~MYPATH/genome_mallotus_dummy.fasta data/genomes/genome_mallotus_dummy.genome
```

And then it can build the database:
```
java -jar snpEff.jar build -gff3 -v genome_mallotus_dummy
```
If you want to, you can look at the database by doing:

```
java -jar ~/Share/ressources/snpEff/snpEff.jar dump genome_mallotus_dummy | less
```
To exit "less" simply press "q"
