#!/bin/bash
###bowtie_mapping.sh 

SAMPLELIST=$1 ### samplelist_all.txt

cd ~/home/physalia_adaptation_course/01_day1/02_genome/
bowtie2-build -f genome_mallotus_dummy.fasta mallotus

mkdir ~/home/mapped/
cd ~/home/mapped/
mkdir samfiles
mkdir bamfiles

for IND in `cat $SAMPLELIST`; do
bowtie2 -p 24 -q --phred33 --end-to-end --very-sensitive --fr --time -x ~/home/physalia_adaptation_course/01_day1/02_genome/mallotus -U ~/home/data/${IND}.fq.gz -S ~/home/mapped/samfiles/${IND}.sam

cd ~/home/mapped/samfiles/
grep -v XS:i: ${IND}.sam > ${IND}.xs.sam
samtools view -b ${IND}.xs.sam | samtools sort - > ../bamfiles/${IND}.bam

done

cd ../
rm -rf samfiles
