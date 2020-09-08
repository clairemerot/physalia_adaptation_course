#!/bin/bash
###stacks_gstacks_2lin.sh
cd
mkdir -p stacks
cd stacks
mkdir -p gstacks_2lin
gstacks -I ~/bamfiles -M ~/scripts/popmap_2lin.txt -O gstacks_2lin -t 3
