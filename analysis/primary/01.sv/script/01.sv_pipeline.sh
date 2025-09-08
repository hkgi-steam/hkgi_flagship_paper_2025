#!/bin/bash

#Configuration
vcf_dir="/mnt/data/patient_data/volume1/genomic_data/latest"
output_dir="../results"
refseq_csv_path="refseq_original.tsv"
sample_id=$1

manta_vcf=$(find "$vcf_dir/$sample_id"/*/ -name "${sample_id}*manta.vcf.gz" | head -n 1)
echo $manta_vcf
cnvkit_vcf=$(find "$vcf_dir/$sample_id"/*/ -name "${sample_id}*cnv.vcf.gz" | head -n 1)
echo $cnvkit_vcf

if [ -z "$manta_vcf" ]; then
    echo "VCF file for patient ID $sample_id not found! Skipping to next sample."
    exit 0
fi

python3 ./sv.py  "$manta_vcf" "$sample_id" "$cnvkit_vcf" "$output_dir"

