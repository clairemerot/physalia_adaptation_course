#!/bin/bash
###stacks_gstacks.sh
mkdir -p stacks
cd stacks
mkdir -p gstacks
cd ~/home/
gstacks -I ~/home/mapped/bamfiles -M physalia_adaptation_course/00_documents/popmap_all.txt -O gstacks -t 24
