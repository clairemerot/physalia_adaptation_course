#!/bin/bash
###bowtie_mapping.sh 

SAMPLELIST=$1 ### samplelist_all.txt

cd ~/Anna/scripts
bowtie2-build -f genome_mallotus_dummy.fasta mallotus

mkdir -p ~/Anna/mapped/
cd ~/Anna/mapped/
mkdir -p samfiles
mkdir -p bamfiles

for IND in `cat $SAMPLELIST`; do
bowtie2 -p 24 -q --phred33 --end-to-end --very-sensitive --fr --time -x ~/Anna/scripts/mallotus -q ~/Anna/data/${IND}.fq.gz > ~/Anna/mapped/samfiles/${IND}.sam >2 bowtie.log

cd ~/Anna/mapped/samfiles/
grep -v XS:i: ${IND}.sam > ${IND}.xs.sam
samtools view -b ${IND}.xs.sam | samtools sort - > ../bamfiles/${IND}.bam

done

cd ../
rm -rf samfiles
