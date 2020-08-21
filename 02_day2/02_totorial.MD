```
vcftools --vcf populations_2lin_random/populations.snps.filter.vcf --site-pi --keep ~/scripts/pop_canada.txt --out pi_canada
vcftools --vcf populations_2lin_random/populations.snps.filter.vcf --site-pi --keep ~/scripts/pop_greenland.txt --out pi_greenland
