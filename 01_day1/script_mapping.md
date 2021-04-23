### Mapping RADseq data to a reference genome
To map our RADseq data to the capelin reference genome we used this script ```physalia_adaptation_course/01_day1/01_scripts/utilities/04_bwa_mem_align_ionproton.sh``` (also on the server at ```~/Share/physalia_adaptation_course/01_day1/01_scripts/utilities/04_bwa_mem_align_ionproton.sh```). 
Here I will provide you with a step-by-step breakdown of that script.
The file starts with the line
```
#!/bin/bash
```
We will write your scripts in bash, and that is the line that tells the server in what language the script is written in.
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
