#!/usr/bin/env bash

set -euo pipefail


if [ "${#@}" -lt 4 ]; then
  echo "Usage: run_eval_syndip.sh <input_dir> <chm_kit_dir> <ref_dir> <ref_name> [filter=true]
Requires Java 8.
Parameters need to be supplied with full path names. 
Input directory needs to have one subdirectory per sample, containing .vcf.gz file(s) for the respective sample. Output summary files will be written here as well."
  exit 1
fi 

input_dir="${1%/}"
chm_kit_dir="${2%/}"
ref_dir="${3%/}"
ref_name="${4}"
filter="${5:-true}"

#ref_name="Homo_sapiens.GRCh38_full_analysis_set_plus_decoy_hla"

export JAVA_HOME=/software/jdk1.8.0_74
export PATH=$JAVA_HOME/bin:$PATH
export PATH=${chm_kit_dir}/:$PATH

if [ ! -f "${ref_dir}/${ref_name}.fa" ]; then
    echo "FASTA file ${ref_dir}/${ref_name}.fa does not exist"
    exit 2
fi

if [ ! -f "${ref_dir}/${ref_name}.fa.fai" ]; then
    echo "FASTA index file ${ref_dir}/${ref_name}.fa.fai does not exist"
    exit 2
fi

SAMPLES=$(ls "${input_dir}")

for sample in ${SAMPLES}; do
    cd "${input_dir}/${sample}"
    for vcf_file in $(ls *.vcf.gz | grep -v flt); do

	if [ ! -d "${ref_dir}/${ref_name}.sdf" ]; then
	    #convert .fa reference to .sdf format required by vcfeval
	    rtg format -o "${ref_dir}/${ref_name}.sdf" "${ref_dir}/${ref_name}.fa"
	fi

	if [[ "${filter}" =~ 'true' ]]; then
	    if [ ! -f "${vcf_file%%.vcf.gz}.flt.vcf.gz" ]; then
		#apply syndip filters
		"${chm_kit_dir}/run-flt" "${vcf_file}"
	    fi

	    #run variant distance based evaluation
	    run-eval -g 38 "${vcf_file%%.vcf.gz}.flt.vcf.gz" | sh

	    #run vcfeval based evaluation for genotype and allele accuracy
	    run-eval -g 38 -s "${ref_dir}/${ref_name}.sdf" "${vcf_file%%.vcf.gz}.flt.vcf.gz" | sh
	else
	    #run variant distance based evaluation    
	    run-eval -g 38 "${vcf_file}" | sh

	    #run vcfeval based evaluation for genotype and allele accuracy    
	    run-eval -g 38 -s "${ref_dir}/${ref_name}.sdf" "${vcf_file}" | sh
	fi
    done
done

