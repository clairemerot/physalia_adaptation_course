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
cp ~/Share/physalia_adaptation_course/01_day1/01_scripts/* ./
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
cd 
mkdir bamfiles
cd bamfiles
ln -s ~/Share/physalia_adaptation_course/01_day1/04_aligned_files/*.bam .
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
###STACKS (http://catchenlab.life.illinois.edu/stacks/) is a software pipeline for building loci from short-read data and for calling variants for population genomics, phylogenmics and linkage maps. STACKS includes two main pipelines: ```denovo_map.pl``` to build a catalog of loci de novo, i.e. it assembles loci from short-read data in absence of a reference genome, and ```ref_map.pl``` to call variants from data mapped to a reference genome. Each of the two pipelines includes several steps.

```denovo_map.pl``` = ```ustacks``` --> ```cstacks``` --> ```sstacks``` --> ```tsv2bam``` --> ```gstacks``` --> ```populations```

```ref_map.pl``` = ```gstacks``` --> ```populations```

As you might have guessed already, the first 4 steps in ```denovo_map.pl``` are to assemble loci from short reads for each individual and to and create a catalog of all loci across the population when you don't have a refernce genome. The next two steps are shared between the two pipeline. ```gstacks``` calls SNPs in each sample and ```populations``` generates population-level summary statistics and input files for a variety of software for downstream analyses.

As most of you have at least a draft reference genome, we'll focus on the ```ref_map.pl``` pipeline. For those interested in the ```denovo_map.pl``` pipeline, this article is a great resource to pick the most appropriate combination of parameters 
```
Paris, J. R., Stevens, J. R., & Catchen, J. M. (2017). 
Lost in parameter space: a road map for stacks. 
Methods in Ecology and Evolution, 8(10), 1360-1373.
```

Although the modules ```gastacks``` and ```populations``` can be run in one go with the ```ref_map.pl``` pipeline, we generally call variants in ```gstacks``` once, and then run ```populations``` several times to tune variant filtering, to print variants for different subsets of samples and to generate F-statistics for each of these subsets.

You should already have all the scripts in your ```~/scripts``` folder.
Let's have a look at the ```stacks_gstacks.sh``` script.
```
#!/bin/bash
###stacks_gstacks.sh
cd
mkdir -p stacks
cd stacks
mkdir -p gstacks
gstacks -I ~/bamfiles -M ~/scripts/popmap_all.txt -O ~/stacks/gstacks -t 3
```
This is simple and self explanatory: we tell the server that we are writing in bash, create a couple of folders (but only if they are not already present with the ```-p``` flag) and run gstacks.

To run this job, we can either use the command ```bash``` before the script or make it executable before running it:
```
nohup bash ~/scripts/stacks_gstacks.sh >& gstacks.log &
```
or
```
chmod +x ~/scripts/stacks_gstacks.sh ### notice that the script name has turned green?
nohup ~/scripts/stacks_gstacks.sh >& gstacks.log &
```

Note that we use ```nohup``` to send the job to the background, so that we have the terminal free for us to do other things and can save the stdout for a file of our choice.

When the job starts you can start checking how it is going from the ```gstacks.log``` file (you can use ```cat```, ```less```, ```more```, etc)

Note 
```
Model: marukilow (var_alpha: 0.01, gt_alpha: 0.05)
```
For most datasets, the authors recommend the default 'marukilow' model, as it takes a Bayesian approach (incorporating information about allele frequencies of the population at each site) and works quite well. The 'marukihigh' model can call more than two alleles per site, but does not use a Bayesian approach. The 'snp' model is the model by Hohenlohe and collaborators and it was the default model in Stacks v1. Below is Maruki + Lynch's paper to fully understand the new models:
```
Maruki, T., & Lynch, M. (2017). 
Genotype Calling from Population-Genomic Sequencing Data. 
G3: Genes|Genomes|Genetics, 7(5), 1393–1404. 
```
To expedite things (and while the first gstacks script runs), we can focus on a smaller dataset. This includes only 80 individuals, 40 from two populations in Canada and 40 from 2 populations in Greenland. Run
```
nohup bash ~/scripts/stacks_gstacks_2lin.sh >& gstacks_2lin.log &
```
Once we have our catalog of variants in the ```gstacks_2lin``` folder we can run populations to apply filters and print variant files in different formats. As ```populations``` is very fast for this small dataset, we can run it from the terminal. But first you make the folder where you want your results to be saved.
```
mkdir populations_2lin_random
populations -P ~/stacks/gstacks_2lin/ -M ~/scripts/popmap_2lin.txt -O populations_2lin_random -t 4 -p 2 -r 0.8 --fstats --vcf --genepop --structure --write-random-snp
```
Let's break the scripts down to understand what ```populations``` is doing here.

-P path to the directory containing the Stacks files.

-M path to a population map

-0 path to a directory where to write the output files

-p minimum number of populations a locus must be present in to process a locus

-r minimum percentage of individuals in a population required to process a locus for that population

--fstats --vcf --genepop --structure for the file formats and stats

--write-random-snp restrict data analysis to one random SNP per locus. This deserves a little more explanation. Most population genetics analyses work under the assumption that your markers are in linkage equilibrium. As you will have multiple SNPs for each locus, and linkage disequilibrium (LD) decays with increasing distance, we select one random SNP for each locus to minimize LD among SNPs. 
It's important that you familiarize with software manuals. Go to <https://catchenlab.life.illinois.edu/stacks/comp/populations.php> and look at the complete list of settings. You can get it also from the terminal if you type ```populations``` (or the program command in general) without additional settings (sometimes ```-h``` or ```--help``` are required to visualize the complete list of settings). 

Also note that the files within each folders have pretty standardized names (check with ```ls -lrth populations_2lin_random```). Pro: we don't need to taylor all our scripts to the individual file names, BUT if you don't specify the right folder you're gonna overwrite the existing files causing a lot of confusion on what setting and samples you used for a particular run.

Although we applied stringent filters to minimize the amount of missing data, these are applied at the SNP level and we could have single individuals that for a variety of reasons (poor DNA quality, low DNA concentration, pipetting errors). We can check the amount of missing data for each individuals with
```
vcftools --vcf populations.snps.vcf --missing-indv --out missing
sort -kr5,5 missing.imiss 
```
In this case % of missing data are below 5%, so we don't need to discard any individual.

#### Congratulations! You have all the input files ready for downstream analyses (plus a bunch of statistics)!

Here you've executed the pipeline on 80 individuals. Before we wrap up today, you need to create the same input files for the full dataset (280 individuals listed in popmap_all.txt) and the Canada dataset (240 individuals listed in popmap_canada.txt). 
First we need to create the catlog of loci from all 280 individuals. Create a new bash script from ```~/scripts/stacks_gstacks_2lin.sh```. What do you need to modify? 
Then, run populations on the new catalog twice, once for all individuals, and once for just the Canadian samples. What do you need to modify here
```
mkdir populations_2lin_random
populations -P ~/stacks/gstacks_2lin/ -M ~/scripts/popmap_2lin.txt -O populations_2lin_random -t 4 -p 2 -r 0.8 --fstats --vcf --genepop --structure --write-random-
```
  





Now, you can run the same analysis on the


Finally, we apply additional filters with the script ```~/scripts/filter_vcf_fast.py```. 
```
chmod +x ~/scripts/filter_vcf_fast.py ### to make it executable
~/scripts/filter_vcf_fast.py ### to check for usage instructions
cd ~/stacks/populations_2lin_random
~/scripts/filter_vcf_fast.py populations.snps.vcf 4 70 0 2 populations.snps.filter.vcf ### you can use this command line for all the different vcf that you generate, as long as you cd to the right folder
```
*** not tested yet
To extract list of SNPs from vcf file for whitelist
cd ~/stacks
bcftools query -f '%CHROM %POS\n' populations_2lin_random/populations.snps.filter.vcf > whitelist_2lin_random.txt

rerun populations to get the updated files with the filetred SNPs
cd ~/stacks
populations -P ~/stacks/gstacks_2lin/ -M ~/scripts/popmap_2lin.txt -O populations_2lin_random_white --fstats --vcf --genepop --structure --whitelist whitelist_2lin_random.txt


populations -P ~/stacks/gstacks/ -M ~/scripts/popmap_all.txt -O populations_all -p 14 -r 0.8 --fstats --vcf --genepop --structure 
populations -P ~/stacks/gstacks_2lin/ -M ~/scripts/popmap_2lin.txt -O populations_2lin -p 2 -r 0.8 --fstats --vcf --genepop --structure 
~/scripts/filter_vcf_fast.py populations.snps.vcf 4 70 0 2 populations.snps.filter.vcf


mkdir populations_canada_random
populations -t 16 -P ~/stacks/gstacks/ -M ~/scripts/popmap_canada.txt -O populations_canada_random -p 12 -r 0.8 --fstats --vcf --genepop --structure --write-random-snp

```
#### to get a list of 1000 random SNPs
grep -v "^#" populations.sumstats.tsv | cut -f 1,4 | sort | uniq | shuf | head -n 1000 | sort -n > whitelist.tsv
grep -v "^#" populations.fst_L-K.tsv | sort -k7,7 | head -n 
cat populations.fst_*.tsv | grep -v "^#" - | sort -rk8,8 | head -n 200 | cut -f 1,6 | sort | uniq | sort -n | head -n 100 > high_fst.whitelist.tsv
mkdir populations_canada_random_white
populations -t 16 -P ~/stacks/gstacks/ -M ~/scripts/popmap_canada.txt -O populations_canada_random_white --fstats --vcf --genepop --structure -W ~/stacks/populations_canada_random/high_fst.whitelist.tsv
