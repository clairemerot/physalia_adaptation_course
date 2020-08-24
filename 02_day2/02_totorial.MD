```
populations -V populations.snps.filter.vcf -M ~/scripts/popmap_2lin.txt -O ../populations_2lin_random_white --fstats --vcf --genepop --structure

mkdir ../populations_all_random_white/
populations -V populations.snps.filter.vcf -M ~/scripts/popmap_all.txt -O ../populations_all_random_white --fstats --vcf --genepop --structure

vcftools --vcf populations_2lin_random/populations.snps.filter.vcf --weir-fst-pop ~/scripts/pop_canada.txt --weir-fst-pop ~/scripts/pop_greenland.txt  --out fst_2lin
vcftools --vcf populations_2lin_random/populations.snps.filter.vcf --weir-fst-pop ~/scripts/pop_canada.txt --weir-fst-pop ~/scripts/pop_greenland.txt  --out fst_2lin_win --fst-window-size 10000 --fst-window-step 10000
vcftools --vcf populations_2lin_random/populations.snps.filter.vcf --site-pi --keep ~/scripts/pop_canada.txt --out pi_canada
vcftools --vcf populations_2lin_random/populations.snps.filter.vcf --site-pi --keep ~/scripts/pop_greenland.txt --out pi_greenland
