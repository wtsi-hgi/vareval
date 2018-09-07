#!/usr/bin/env bash                                                                                                                                                                                                                                                                                                           

set -euo pipefail

if [ "${#@}" -lt 3 ]; then
  echo "Usage: split_truth_by_af.sh <input_vcf> <AC_strat_file> <AF_strat_file>
Parameters need to be supplied with full path names.
<AC_strat_file> needs to contain desired allele counts, 1 per line.
<AF_strat_file> needs to contain the boundaries for desired allele frequency bins, 1 per line."
  exit 1
fi

input_vcf="${1}"
ac_file="${2}"
af_file="${3}"
region="${4:-none}"


if [ ! -f "${input_vcf}" ]; then
    echo "Input vcf file does not exist."
    exit 2
fi

input_vcf_name=$(basename "${input_vcf}")

if [[ "${region}" =~ 'none' ]]; then
    region_flag=""
    region="wg"
else
    region_flag=" -r ${region}"
fi

if [ ! -f "${ac_file}" ]; then
    echo "AC stratification file does not exist."
    exit 2
else
    for ac in $(cat "${ac_file}"); do
	bcftools view -Oz"${region_flag}" -o "${input_vcf_name%.vcf.gz}.${region}.AC${ac}.vcf.gz" -e AC!="${ac}" "${input_vcf}"
    done 
fi

if [ ! -f "${af_file}" ]; then
    echo "AF stratification file does not exist."
    exit 2
else
    af_list=( $(cat "${af_file}") ) 
    for ind in "${!af_list[@]}"; do 
	if [[ "${ind}" == 0 ]]; then
	    if [[ "${af_list[${ind}]}" != 0 ]]; then
		bcftools view -Oz"${region_flag}" -o "${input_vcf_name%.vcf.gz}.${region}.AF0-${af_list[${ind}]}.vcf.gz" -e "AF>=${af_list[${ind}]}" "${input_vcf}"
	    fi
	else
            bcftools view -Oz"${region_flag}" -o "${input_vcf_name%.vcf.gz}.${region}.AF${af_list[${ind}-1]}-${af_list[${ind}]}.vcf.gz" -e "AF>=${af_list[${ind}]} || AF<${af_list[${ind}-1]}"  "${input_vcf}"
	    if [[ "${ind}" == "$((${#af_list[@]}-1))" ]] && [[ "${af_list[${ind}]}" != 1 ]]; then
		bcftools view -Oz"${region_flag}" -o "${input_vcf_name%.vcf.gz}.${region}.AF${af_list[${ind}]}-1.vcf.gz" -e "AF<${af_list[${ind}]}" "${input_vcf}"
	    fi
        fi
    done 
fi