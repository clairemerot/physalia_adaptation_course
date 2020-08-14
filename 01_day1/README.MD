# Day1: Handling NGS data: From raw reads to SNPs matrix

## Step 1 : XX

### dacdca

Welcome to the practical portion of the Adaptation Genomics course!

Here are the instructions to log in the AWS server.

Once you have successfully logged in, the first thing you will do is to create your own folder and generate a symbolic link so that we don't need to personalize all the scripts
```
mkdir *yourlastname*
ln -s ~/*yourlastname* ~/home
```

Copy this repository in your folder
```
cd ~/home
git clone  https://github.com/clairemerot/physalia_adaptation_course.git
```


## Step 1 : Mapping short reads to the capelin (dummy) reference genome
The data we'll be working on are from capelin (Mallotus villosus), a small marine fish, from Cayuela et al. 2020, Molecular Ecology. They were generated on the IonTorrent platform but all the preprocessing (demultiplexing, adapter removal, read trimming, quality filtering) has already been done. You can find example script for these preliminary steps here and here.

To run analysis that take some time we will submit jobs through bash scripts.
To map the filtered data to the capelin reference genome, make the script ```bowtie_mapping.sh``` in folder xxx executable and submit it with
```
chmod +x physalia_adaptation_course/01_day1/01_scripts/bowtie_mapping.sh
physalia_adaptation_course/01_day1/01_scripts/bowtie_mapping.sh physalia_adaptation_course/00_documents/samplelist_all.txt
```
To run a bash script like this you can also type
```
bash physalia_adaptation_course/01_day1/01_scripts/bowtie_mapping.sh physalia_adaptation_course/00_documents/samplelist_all.txt
``` 
without making it executable.

### File formats
While we wait we can explore the read file (.fastq format) and the alignment file (.sam/.bam format), which are the input and output of the alignemnt, respectively.
#### FASTQ
```
@NB551191:35:HNWGCBGX3:1:11101:9622:1037 1:N:0:AGGCAGAA+ACTCTAGG
CCTTGNTGCACCTGTGGCATGGAGACCGAATCTTGTGGGGAAACAATCATTTCTTCAGGTCTGAGCTCTCAGATTT
+
AAAAA#EAEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEAE/E/AEEEEEEEEEE/EEEEEEEEEEEEEEEEEEE
@NB551191:35:HNWGCBGX3:1:11101:20086:1038 1:N:0:AGGCAGAA+ACTCTAGG
CTGGTNCACCATCCTTGTGTGCTGTTTCATGACAGTAATTACTGAGAGGGTCTGCAATTCAGATCACCTGAAACTC
+
AAAAA#EEEEAEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
@NB551191:35:HNWGCBGX3:1:11101:4269:1051 1:N:0:AGGCAGAA+ACTGTAGG
TAACAGCACAGAGGATTGAATAAGGTGAGAGCAAAAGTCCTACTACTTATTCAGGCCCCATGTAGCAGTATTCCTC
+
AAAAA6EEEEEEEEEEEEEA/EEEEEEEAEEEEEEE6/EA/EEEEEE/E//EA//EEEEEEEE6A/E6/EE/EEEE
```
The snippet above shows 3 reads, with each read coming with 4 lines of information.
1. Sequence ID and information
2. The actual read sequence
3. Generally plus a '+' sign but it could contain additional information
4. The quality scores for each base in line 2

Knowing this, we can calculate the number of reads in a file by counting the number of lines
```
wc -l file.fastq
```
and divide by 4.

When we have paired reads, they are generally stored in two separate files and are paired by position: the first read in ```read1.fastq``` is paired with the first read in ```read2.fastq```, the second read in ```read1.fastq``` is paired with the second read in ```read2.fastq```, and so on. One way to check if your files are intact is to ensure that the two read files have the same amount of lines/reads.
#### SAM/BAM
The alignment file contains much more information. In addition to all the information contained in the .fastq files, we now have information of the quality of the alignment and the position of where those reads mapped in the genome. Although we produce .sam files with the alignment, we quickly convert them to .bam, which is in binary format. These files are not readable directly (but you can use ```samtools view```).

Back to mapping...

