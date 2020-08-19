# Day1: Handling NGS data: From raw reads to SNPs matrix
## Step 1: Getting familiar with the UNIX environment
Although you should have familiarized with working on a UNIX environment before the beginning of the course, here below are the commands that we will use most commonly. After you have logged in, create your own folder and move in it. You will run all your scripts from here
```
mkdir yourname
cd yourname
```
Then copy some of the files we will need today in a new folder
```
mkdir scripts
cd scripts
cp ~/Share/physalia_adaptation_course/00_documents/* ./ ###this command copies all the files in the folder 00_documents into the folder you're currently in
```
Use ```ls``` to get a list of the files we have copied over, and look at the content of one of those files with
```
less popmap_2lin.txt 
```
You can stroll down by pressing the bar on your keyboard and ```q``` to exit from the screen.

If you want to print a few lines to the screen (which won't disapper) type
```
head popmap_2lin.txt
```
which will print the first 10 lines by default. You can specify the number of lines with ```-n 25```. Same for ```tail```, which shows the end of the file.

To edit scripts and files, we'll use ```nano``` as text editor
Type 
```
nano popmap_2lin.txt
``` 
and move around using the arrows on your keyboard.
At the bottom of your screen are the shortcuts for different commands. For example, to modify the file and save it with a different name, write something in the first line, press ctlr + x and edit the name of the file into ```test_nano.txt```. Hit enter and the file will be saved.
Now, let's remove this file with
```
rm test_nano.txt
```
WARNING!!! ```rm``` removes files and folders for good, and they can't be retrieved. So please keep that in mind!

## Step 2: Getting from raw data to mapped data
### Pre-processing of raw data
The data we'll be working on are from capelin (Mallotus villosus), a small marine fish, from Cayuela et al. 2020, Molecular Ecology. They were generated on the IonTorrent platform but all the preprocessing (demultiplexing, adapter removal, read trimming, quality filtering) has already been done. You can find example script for these preliminary steps in the scripts folder of day 1.

The data, raw or cleaned, are stored in .fastq files and here is what they look like.

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

### Mapping data to a reference genome
As we're analyzing a large dataset that would take many hours to align to the reference genome, we also did the mapping for you using this script ```physalia_adaptation_course/01_day1/01_scripts/utilities/04_bwa_mem_align_ionproton.sh```. 
Here I will provide you with a step-by-step breakdown of that script.
The file starts with the line
```
#!/bin/bash
```
We will write your scripts in bash, and that is the line that tell the server in what language the script is written in.
A few lines starting with ```#SBATCH``` follow. We won't need these lines because we won't use a job scheduler, but this is a good example for the people that have SLURM on their institution's servers.

We set a few variables as shown below, These are very handy to make this script easy to follow and more customizable. For example, the number of CPUs is not hard-coded in the script and can be conveniently changed when submitting the script.
```
# Global variables
GENOMEFOLDER="02_genome"
GENOME="genome_mallotus_dummy.fasta"
DATAFOLDER="03_raw_reads"
ALIGNEDFOLDER="04_aligned_files"
NCPU=$1
```

First, it uses a for loop in bash to loop through the list of fastq files to map to the reference genome
```
for file in $(ls -1 "$DATAFOLDER"/*.fq.gz)
do
```
then, it uses a one-liner (even though it looks like multiple lines for the use of backslashes \) to map the data and pipe (with |) the .sam output into the binary format .bam
```
bwa mem -t "$NCPU" -k 19 -c 500 -O 0,0 -E 2,2 -T 0 \
        -R "$ID" \
        "$GENOMEFOLDER"/"$GENOME" "$DATAFOLDER"/"$name" 2> /dev/null |
        samtools view -Sb -q 1 -F 4 -F 256 -F 2048 \
        - > "$DATAFOLDER"/"${name%.fq.gz}".bam
```
The alignments in bam files are sorted and indexed
```
samtools sort --threads "$NCPU" -o "$DATAFOLDER"/"${name%.fq.gz}".sorted.bam \
        "$DATAFOLDER"/"${name%.fq.gz}".bam

    samtools index "$DATAFOLDER"/"${name%.fq.gz}".sorted.bam
```
and the intermediate files are deleted with
```
rm "$DATAFOLDER"/"${name%.fq.gz}".bam
```

Now, move to your home folder, make a new folder for the alignment files and create a symbolic link to the bamfiles stored in the```~/Share``` folder
```
cd ~/Anna
mkdir bamfiles
cd bamfile
ln -l ~/Share/bamfiles/* ./
```

#### SAM/BAM
The alignment file contains much more information than the raw data .fastq files. In addition to all the information contained in the .fastq files, we now have information of the quality of the alignment and the position of where those reads mapped in the genome. Although we produce .sam files with the alignment, we quickly convert them to .bam, which is in binary format. These files are not readable directly but you can use 
```
samtools view -h xxx.bam | less 
```
With the pipe sign ```|```, we pipe the output of the first command into the second command.
The first few lines starting with ```@``` is the header. If you have a lot of scaffolds in your reference the header can be very long. If you want to see the first few alignments remove ```-h``` from the command above.

nohup ./model >outfile.txt &



        
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

## Step 3: Call variants in STACKS
STACKS (http://catchenlab.life.illinois.edu/stacks/) is a software pipeline for building loci from short-read data and for calling variants for population genomics, phylogenmics and linkage maps. STACKS includes two main pipelines: ```denovo_map.pl``` to build a catalog of loci de novo, i.e. it assembles loci from short-read data in absence of a reference genome, and ```ref_map.pl``` to call variants from data mapped to a reference genome. Each of the two pipelines includes several steps.

```denovo_map.pl``` = ```ustacks``` --> ```cstacks``` --> ```sstacks``` --> ```tsv2bam``` --> ```gstacks``` --> ```populations```

```ref_map.pl``` = ```gstacks``` --> ```populations```

As you might have guessed already, the first 4 steps in ```denovo_map.pl``` are to assemble loci from short reads for each individual and to and create a catalog of all loci across the population. The next two steps are shared between the two pipeline. ```gstacks``` call SNPs in each sample and ```populations``` generates population-level summary statistics and input files for a variety of software for downstream analyses.
