#!/bin/bash
###stacks_gstacks.sh
cd
mkdir -p stacks
cd stacks
mkdir -p gstacks
gstacks -I ~/bamfiles -M ~/scripts/popmap_all.txt -O gstacks -t 3
