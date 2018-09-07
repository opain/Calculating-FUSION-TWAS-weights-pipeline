#!/bin/bash

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    --PLINK_prefix)
    target_plink="$2"
    shift # past argument
    shift # past value
    ;;
    --phenotype_file)
    phenotype_file="$2"
    shift # past argument
    shift # past value
    ;;
    --coordinate_file)
    coordinate_file="$2"
    shift # past argument
    shift # past value
    ;;
    --plink)
    plink="$2"
    shift # past argument
    shift # past value
    ;;
    --gcta)
    gcta="$2"
    shift # past argument
    shift # past value
    ;;
    --gemma)
    gemma="$2"
    shift # past argument
    shift # past value
    ;;
    --ld_ref_dir)
    ld_ref_dir="$2"
    shift # past argument
    shift # past value
    ;;
    --fusion_software)
    fusion_software="$2"
    shift # past argument
    shift # past value
    ;;
    --output_dir)
    output="$2"
    shift # past argument
    shift # past value
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

echo target_plink = "${target_plink}"
echo phenotype_file = "${phenotype_file}"
echo coordinate_file = "${coordinate_file}"
echo plink = "${plink}"
echo gcta = "${gcta}"
echo gemma = "${gemma}"
echo ld_ref_dir = "${ld_ref_dir}"
echo fusion_software = "${fusion_software}"
echo output = "${output}"

feature_name=$(awk "NR==${SGE_TASK_ID}" <(cut -f 4 -d ' ' ${coordinate_file}))

echo feature_name = "${feature_name}"

Rscript ./OP_TWAS_weights_using_fusion.R \
--PLINK_prefix ${target_plink} \
--phenotype_file ${phenotype_file} \
--coordinate_file ${coordinate_file} \
--gene_name ${feature_name} \
--plink ${plink} \
--gcta ${gcta} \
--gemma ${gemma} \
--ld_ref_dir ${ld_ref_dir} \
--fusion_software ${fusion_software} \
--output_dir ${output}
