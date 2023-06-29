# Day 1 - tutorial 2: Handling NGS data: From raw reads to SNPs matrix

## Step 1: Getting familiar with the UNIX environment
Although you should have familiarized with working on a UNIX environment before the beginning of the course, here below are the commands that we will use most commonly. When you log in, you're in your home folder. You can go back here any time by typing ```cd``` (and pressing return).

First, we'll create a folder and copy some of the files we will need today in this new folder
```
mkdir scripts
cd scripts
cp ~/Share/physalia_adaptation_course/00_documents/* . ###this command copies all files (*) that were in the folder 00_documents into the folder you're currently in
cp -r ~/Share/physalia_adaptation_course/01_day1/01_scripts . ###this command copies a folder (note the -r juste after the command "cp" in the folder 00_documents into the folder you're currently in
```
Use ```ls``` to get a list of the files we have copied over, and look at the content of one of those files with the ```less``` command.
```
less 00_documents/popmap_2lin.txt 
```
Or another possibility is to navigate inside the folders with ```cd WantedFolder``` command and then use ```less WantedFile``` to display the file content.

Example for the population map : ```cd 00_documents/popmap_2lin.txt```

Then, you can stroll down by pressing the bar on your keyboard and exit from the screen with ```q```.

If you want to print a few lines to the screen (which won't disapper) type
```
head 00_documents/popmap_2lin.txt
```
which will print the first 10 lines by default. You can specify the number of lines with ```-n 25```. Same for ```tail```, which shows the end of the file.

To edit scripts and files, we'll use ```nano``` as text editor.
Type 
```
nano popmap_2lin.txt
``` 
and move around using the arrows on your keyboard.
At the bottom of your screen are the shortcuts for different commands. For example, to modify the file and save it with a different name, write something in the first line, press ```ctlr + x``` and edit the name of the file into ```test_nano.txt```. Hit enter and the file will be saved.
Now, let's remove this file with
```
rm test_nano.txt
```
WARNING!!! ```rm``` removes files and folders for good, and they can't be retrieved. So please keep that in mind!

## Step 2: Getting from raw data to mapped data
### Pre-processing of raw data
The data we'll be working on are from capelin (*Mallotus villosus*), a small marine fish, from Cayuela et al. 2020, Molecular Ecology. They were generated on the IonTorrent platform but all the preprocessing (demultiplexing, adapter removal, read trimming, quality filtering) has already been done. You can find example script for these preliminary steps in the scripts folder of day 1.

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
zcat -f  ~/Share/physalia_adaptation_course/01_day1/03_raw_reads/A_01.fq.gz | wc -l 
```
and divide by 4.
Or 
```
zcat -f  ~/Share/physalia_adaptation_course/01_day1/03_raw_reads/A_01.fq.gz | sed -n '2~4p' | wc -l
```
directly from compressed fastq files

When we have paired reads, they are generally stored in two separate files and are paired by position: the first read in ```read1.fastq``` is paired with the first read in ```read2.fastq```, the second read in ```read1.fastq``` is paired with the second read in ```read2.fastq```, and so on. One way to check if your files are intact is to ensure that the two read files have the same amount of lines/reads.

### Mapping data to a reference genome
As we're analyzing a large dataset that would take many hours to align to the reference genome, we also did the mapping for you. 
To find the script that we used and a step-by-step breakdown of that script, click on the link below.

[Script to map to a reference genome](./script_mapping.md)

Now, move to your home folder, make a new folder for the alignment files and create a symbolic link to the bamfiles stored in the```~/Share``` folder. With a symbolic link, you add the files to a path of your choice, but because they are not physically there you're not taking up the space.
```
cd 
mkdir bamfiles
cd bamfiles
ln -s ~/Share/bamfiles/*.bam .
```

#### SAM/BAM
The alignment file contains much more information than the raw data .fastq files. In addition to all the information contained in the .fastq files, we now have information of the quality of the alignment and the position of where those reads mapped in the genome. Although we produce .sam files with the alignment, we quickly convert them to .bam, which is in binary format. These files are not readable directly but you can use 
```
samtools view -h A_01.bam | less 
```
With the pipe sign ```|```, we pipe the output of the first command into the second command.
The first few lines starting with ```@``` is the header. If you have a lot of scaffolds in your reference the header can be very long. If you want to see the first few alignments remove ```-h``` from the command above.

(exit less by pressing q)
        
## Step 3: Call variants in STACKS from RADseq data
STACKS <http://catchenlab.life.illinois.edu/stacks/> is a bioinformatic pipeline for building loci from short-read data from GBS/RADseq libraries and calling variants for population genomics, phylogenomics and linkage maps. STACKS includes two main pipelines: ```denovo_map.pl``` to build a catalog of loci *de novo*, i.e. it assembles loci from short-read data in absence of a reference genome, and ```ref_map.pl``` to call variants from data mapped to a reference genome. Each of the two pipelines includes several steps. Note that everything here is based on STACKS v.2.5.3.

```denovo_map.pl``` = ```ustacks``` --> ```cstacks``` --> ```sstacks``` --> ```tsv2bam``` --> ```gstacks``` --> ```populations```

```ref_map.pl``` = ```gstacks``` --> ```populations```

As you might have guessed already, the first 4 steps in ```denovo_map.pl``` are for assembling loci from short reads for each individual and for creating a catalog of all loci across the population when you don't have a reference genome. The next two steps are shared between the two pipeline. ```gstacks``` builds a catalog of loci and calls SNPs in each sample and ```populations``` generates population-level summary statistics and input files for a variety of software for downstream analyses.

As most of you have at least a draft reference genome, we'll focus on the ```ref_map.pl``` pipeline. For those interested in the ```denovo_map.pl``` pipeline, this article is a great resource to pick the most appropriate combination of parameters 
```
Paris, J. R., Stevens, J. R., & Catchen, J. M. (2017). 
Lost in parameter space: a road map for stacks. 
Methods in Ecology and Evolution, 8(10), 1360-1373.
```

Although the modules ```gstacks``` and ```populations``` can be run in one go with the ```ref_map.pl``` pipeline, we generally call variants in ```gstacks``` once, and then run ```populations``` several times to tune variant filtering, to print variants for different subsets of samples and to generate F-statistics for each of these subsets.

Please copy the scripts inside your folder and make executable the one we are interested in 
```
cd # come back to your home
cd scripts #get inside your scripts
#now copy files inside your 01_scripts folder. (replace the one you may have copy earlier as I edited them)
cp ~/Share/physalia_adaptation_course/01_day1/01_scripts/*.sh 01_scripts/
#make the file executable
chmod +x ~/scripts/01_scripts/stacks_gstacks.sh ### notice that the script name has turned green?
#open the file
less ~/scripts/01_scripts/stacks_gstacks.sh
#to exit, press 'q'
```

Let's have a look at the ```stacks_gstacks.sh``` script.
```
#!/bin/bash
###stacks_gstacks.sh
cd
mkdir -p stacks
cd stacks
mkdir -p gstacks
gstacks -I ~/bamfiles -M ~/scripts/00_documents/popmap_all.txt -O ~/stacks/gstacks -t 3
```
This is simple and self explanatory: we tell the server that we are writing in bash, create a couple of folders (but only if they are not already present with the ```-p``` flag) and run ```gstacks```.

To run this job, we can either use the command ```bash``` before the script or make it executable before running it:
```
nohup bash ~/scripts/01_scripts/stacks_gstacks.sh >& gstacks.log &
```
or
```

nohup ~/scripts/01_scripts/stacks_gstacks.sh >& gstacks.log &
```

Note that we use ```nohup``` to send the job to the background, so that we have the terminal free for us to do other things and can save the stdout to a file of our choice.

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
nohup bash ~/scripts/01_scripts/stacks_gstacks_2lin.sh >& gstacks_2lin.log &
```
Once we have our catalog of variants in the ```gstacks_2lin``` folder we can run ```populations``` to apply filters and print variant files in different formats. As ```populations``` is very fast for this small dataset, we can run it from the terminal. But first you make the folder where you want your results to be saved.
```
cd ~/stacks
mkdir populations_2lin_random
populations -P ~/stacks/gstacks_2lin/ -M ~/scripts/00_documents/popmap_2lin.txt -O ~/stacks/populations_2lin_random -t 4 -p 2 -r 0.8 --fstats --vcf --genepop --structure --write-random-snp
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

Although we applied stringent filters to minimize the amount of missing data, these are applied at the SNP level and we could have single individuals that for a variety of reasons (e.g., poor DNA quality, low DNA concentration, pipetting errors) have a high proportion of missing data. We can check the amount of missing data for each individuals with
```
cd populations_2lin_random
vcftools --vcf populations.snps.vcf --missing-indv --out missing
sort -rk5,5 missing.imiss | less # with this you will visualize individuals by % of missing data, from the highest to the lowest
```
In this case % of missing data are below 5%, so we don't need to discard any individual.

#### Congratulations! You have all the input files ready for downstream analyses (plus a bunch of statistics) for the reduced dataset!

As you can see if you  do ```ls``` to list the files in one of the folders storing the output from your ```populations``` runs, STACKS produces many output files. The manual provides a thorough explanation of all the different file formats. Go to <https://catchenlab.life.illinois.edu/stacks/manual/> and scroll down to section 6.6. For each output type, open the corresponding file on the server, and take your time to explore them.

Here you've executed the pipeline on 80 individuals. Before we wrap up today, you need to create the same input files for the full dataset (280 individuals listed in popmap_all.txt) and the Canada dataset (240 individuals listed in popmap_canada.txt). 
You should have already run the script ```~/scripts/gstacks_stacks.sh```. 
Now, run ```populations``` on the new catalog twice, once for all individuals, and once for just the Canadian samples. What do you need to modify here?
```
cd ~/stacks
mkdir populations_2lin_random
populations -P ~/stacks/gstacks_2lin/ -M ~/scripts/00_documents/popmap_2lin.txt -O ~/stacks/populations_2lin_random -t 4 -p 2 -r 0.8 --fstats --vcf --genepop --structure --write-random-snp
```
 
Make sure you print these outputs in each of these two folders ```populations_all_random``` and ```populations_canada_random```.

If you are struggling with this last step, have a look at the solution [here](./populations_script.md).

## Step 4. Variant calling from whole genome resequencing data
First, let's set up new folders for these analyses
```
cd
mkdir wgr #for whole genome resequencing
mkdir wgr/snps_bcftools
```

### SNPs calling
Copy bam files (symbolic links) and their indexes storing aligned reads to reference genome
```
cd wgr
ln -s /home/ubuntu/Share/WGS_bam/*bam* .
ls *.bam > bams_list.txt ### create list of bam files
```
Call variants using bcftools
```
bcftools mpileup -Ou -f ~/Share/ressources/genome_mallotus_dummy.fasta -b bams_list.txt -r Chr1 -q 5 -I -a AD,DP,SP,ADF,ADR -d 200 | bcftools call - -mv -Ov > snps_bcftools/capelin_wgs_unfiltered.vcf

```
Variant filtering is done in two steps. First we use 'bcftools filter' to filter SNPs based on mapping quality. Then, we apply a series of filters using VCFtools. Note that I'm calling variants from chromosome 1 only to save time. 
```
cd snps_bcftools
bcftools filter -e 'MQ < 30' capelin_wgs_unfiltered.vcf -Ov > capelin_wgs_filtered.tmp.vcf
### Count number of unfiltered SNPs surviving filters
grep -v ^\#\# capelin_wgs_filtered.tmp.vcf | wc -l  

vcftools --vcf capelin_wgs_filtered.tmp.vcf \
    --max-alleles 2 \
    --max-missing 0.7 \
    --minDP 3 \
    --maf 0.05 \
    --minQ 30 \
    --minGQ 20 \
    --recode \
    --stdout > capelin_wgs_filtered.vcf 


### Count number of SNPs surviving filters
grep -v "#" capelin_wgs_filtered.vcf | wc -l 
```

Because the variant calling takes quite some time, you can copy the original unfiltered and clean files to your working directory to explore the different files and play with filtering.
```
cp ~/Share/WGS_bam/snps_bcftools/capelin_wgs_*.vcf ~/wgr/snps_bcftools/.
```
What is the proportion of SNPs surviving filters? How many are left? According to your dataset and needs, you may want to change or tweak these filters.

Finally, prepare the data for downstream analyses by recoding the genotypes in 012 format with 
```
vcftools --vcf capelin_wgs_filtered.vcf --012 --out capelin_wgs_filtered
```
and download the resulting files on your computer.

### The end.
