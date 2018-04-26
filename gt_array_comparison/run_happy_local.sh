# Script to run hap.py locally

set -euo pipefail

#cd tests

bgzip -c truth.vcf > truth.vcf.gz
bgzip -c calls.vcf > calls.vcf.gz

tabix -f -p vcf truth.vcf.gz
tabix -f -p vcf calls.vcf.gz

cat truth.vcf | grep -v "#" | awk 'BEGIN { OFS="\t"; } { print $1, $2-1, $2-1+length($4); }' - > truth.bed

hap.py truth.vcf.gz calls.vcf.gz -r ../reference.fa --write-vcf -o out/out --engine=vcfeval -f truth.bed

zgrep -v "UNK" out/out.vcf.gz > out/out_no_unk.vcf