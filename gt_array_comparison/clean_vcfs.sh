#!/usr/bin/env bash

set -o pipefail

set -x
'''
module purge
module add hgi/tabix/git-1ae158a
module add hgi/vcftools/0.1.14
module add hgi/bcftools/1.6-htslib-1.6-htslib-plugins-6f2229e0-irods-git-4.2.2-plugin_kerberos-2.0.0
module add hgi/parallel/20140922
'''
USAGE=$(cat -<< EOF
USAGE: $0 --reference/-r <FILE> --sample_list/-s <FILE> [--exome_mask/-e <FILE>] <VCF_FILE>
EOF
)

reference=""
sample_list=""
exome_mask=""
input_vcf=""

OPTIONS=r:s:e:i:
LONGOPTIONS=reference:,sample_list:,exome_mask:,input_vcf:

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
        -i|--input_vcf)
            input_vcf="$2"
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
            echo $1
            exit 3
            ;;
    esac
done

if [[ -z "$input_vcf" ]]; then
    echo "$0: The parameter --input_vcf is required."
    exit 4
fi

if [[ -z "$reference" ]]; then
    echo "$0. The parameter --reference is required."
    exit 4
fi

if [[ -z "$sample_list" ]]; then
    echo "$0. The parameter --sample_list is required."
    exit 4
fi


# # default is ELGH_ROH.gvcf_to_vcf.20180206.hard-filtered.reheaded.vcf.gz
# input_file=$1
# # default is /lustre/scratch118/humgen/hgi/users/mercury/2017-2018-variant-caller-eval/dragen/ref-v6/Homo_sapiens.GRCh38_full_analysis_set_plus_decoy_hla.fa
# fasta_reference=$2
# # default is ../exon_regions_chr
# exome_mask=$3
# # default is ./111_wes_sample_list
# sample_list=$4

alias my_parallel="parallel --no-notice"

input_file_base=$(basename ${input_vcf} .vcf.gz)
working_dir=$(dirname ${input_vcf})

# index input vcf
tabix -p vcf ${input_vcf}

# subset columns
bcftools view -Oz -S ${sample_list} ${input_vcf} > ${input_file_base}.subset.vcf.gz

# index subset vcf
tabix -p vcf ${input_file_base}.subset.vcf.gz

# split vcf by chrom
parallel "tabix -h ${input_file_base}.subset.vcf.gz chr{} | bgzip -c > ${input_file_base}.chr{}.vcf.gz" ::: {1..22}

# index chrom-level vcfs
parallel "bcftools index -t ${input_file_base}.chr{}.vcf.gz" ::: {1..22}

if [[ -z "$exome_mask" ]]; then
    exome_mask_option=" -R ${exome_mask}"
else
    exome_mask_option=""
fi

# normalize vcf + apply exome mask (-R ../exon_regions_chr flag; should be optional)
parallel "bcftools norm -c s -f ${reference} -m +any ${input_file_base}.chr{}.vcf.gz | bcftools annotate -Oz -x FORMAT/AD > ${input_file_base}.chr{}.normed.noAD.vcf.gz" ::: {1..22}

rm *chr{1..22}.vcf.gz*
# we should not remove the input file
#rm ELGH_ROH.gvcf_to_vcf.20180206.hard-filtered.vcf.gz*

#ls ./${input_file_base}.chr*.vcf.gz*

# more intermediate file removal
#(rm ./${input_file_base}.chr*.vcf.gz* || true)

# index vcf obtained in prev step
parallel "bcftools index -t ${input_file_base}.chr{}.normed.noAD.vcf.gz" ::: {1..22}

