#!/usr/bin/env bash

set -euf -o pipefail

set -x

module purge
module add hgi/tabix/git-1ae158a
module add hgi/vcftools/0.1.14
module add hgi/bcftools/1.6-htslib-1.6-htslib-plugins-6f2229e0-irods-git-4.2.2-plugin_kerberos-2.0.0
module add hgi/parallel/20140922

alias my_parallel="parallel --no-notice"

USAGE=$(cat -<< EOF
USAGE: $0 --reference/-r <FILE> --sample_list/-s <FILE> [--exome_mask/-e <FILE>] <VCF_FILE>
EOF
)

reference=""
sample_list=""
exome_mask=""

OPTIONS=r:s:e:
LONGOPTIONS=reference:,sample_list:,exome_mask:

PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTIONS --name "$0" -- "$@")
if [[ $? -ne 0 ]]; then
    exit 2
fi
eval set -- "$PARSED"

while true; do
    case "$1" in
        -r|--reference)
            reference="$2"
            shift 2
            ;;
        -s|--sample_list)
            sample_list="$2"
            shift 2
            ;;
        -e|--exome_mask)
            exome_mask="$2"
            shift 2
            ;;
        -h|--help)
            echo ${USAGE}
            exit 0
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "Programming error"
            exit 3
            ;;
    esac
done

if [[ $# -ne 1 ]]; then
    echo "$0: An input VCF is required."
    exit 4
fi
input_vcf=$1

if [[ -z "$reference" ]]; then
    echo "$0. The parameter --reference needs to be defined."
    exit 4
fi

if [[-z "$sample_list"]]; then
    echo "$0. The parameter --sample_list needs to be defined."
    exit 4
fi

input_file_base=$(basename ${input_file} .vcf.gz)

# index input vcf
tabix -p vcf ${input_file}

# subset columns
bcftools view -Oz -S ${sample_list} ${input_file_base}.vcf.gz > ${input_file_base}.subset.vcf.gz

# index subset vcf
tabix -p vcf ${input_file_base}.subset.vcf.gz

# split vcf by chrom
my_parallel "tabix -h ${input_file_base}.subset.vcf.gz chr{} | bgzip -c > ${input_file_base}.chr{}.vcf.gz" ::: {1..22}

# index chrom-level vcfs
my_parallel "tabix -p vcf ${input_file_base}.chr{}.vcf.gz" ::: {1..22}

if [[ -z "$exome_mask" ]]; then
    exome_mask_option=" -R ${exome_mask}"
else
    exome_mask_option=""
fi


# normalize vcf + apply exome mask (-R ../exon_regions_chr flag; should be optional)
my_parallel "bcftools norm -c s -Oz -f ${reference} -m -any${exome_mask_option} ${input_file_base}.chr{}.vcf.gz > ${input_file_base}.chr{}.exome.normed.vcf.gz" ::: {1..22}

# we should not remove the input file
#rm ELGH_ROH.gvcf_to_vcf.20180206.hard-filtered.vcf.gz*

# more intermediate file removal
rm ${input_file_base}.chr*.vcf.gz*

# retain only PASS columns
my_parallel "bcftools view -Oz -f PASS ${input_file_base}.chr{}.exome.normed.vcf.gz > ${input_file_base}.chr{}.exome.normed.only-pass.vcf.gz" ::: {1..22}

# index vcf obtained in prev step
my_parallel "tabix -p vcf ${input_file_base}.chr{}.exome.normed.only-pass.vcf.gz" ::: {1..22}
