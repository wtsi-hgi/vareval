# Script to run hap.py locally

set -euo pipefail

cd tests

bgzip -c truth.vcf > truth.vcf.gz
bgzip -c calls.vcf > calls.vcf.gz

tabix -f -p vcf truth.vcf.gz
tabix -f -p vcf calls.vcf.gz

hap.py truth.vcf.gz calls.vcf.gz -r reference.fa --write-vcf -o out/out --engine=vcfeval

bgzip -d -f out/out.vcf.gz