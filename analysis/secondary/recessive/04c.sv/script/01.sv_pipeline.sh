#!/bin/bash
#Configuration
vcf_dir="/mnt/data/patient_data/volume1/genomic_data/latest"
projectDir="/home/_shared/jscliu/project/2025/Flagship/analysis/secondary/recessive/04c.sv"
output_dir="${projectDir}/results"
refseq_csv_path="${projectDir}/script/refseq_original.tsv"
bed_file="/home/_shared/jscliu/project/2025/Flagship/reference/hkgi_recessive.ext200bp.bed"

# Input
sample_id=$1
manta_vcf=$(find "$vcf_dir/$sample_id"/*/ -name "${sample_id}*.manta.vcf.gz" | head -n 1)
echo $manta_vcf
cnvkit_vcf=$(find "$vcf_dir/$sample_id"/*/ -name "${sample_id}*.cnv.vcf.gz" | head -n 1)
echo $cnvkit_vcf

if ls "${output_dir}/${sample_id}"*sv.tsv 1> /dev/null 2>&1; then
    echo "SV file for patient ID ${sample_id} already exists. Skipping sample."
    exit 0
fi

if [ -z "$manta_vcf" ]; then
    echo "VCF file for patient ID $sample_id not found! Skipping to next sample."
    exit 0
fi

subset_manta_vcf="$output_dir/${sample_id}_carrier_manta.vcf.gz"
subset_cnvkit_vcf="$output_dir/${sample_id}_carrier_cnv.vcf.gz"

apptainer run /mnt/bioinformatics_tools/tools/bcftools-1.18.sif bcftools view -R "$bed_file" "$manta_vcf" -Oz -o "$subset_manta_vcf"
apptainer run /mnt/bioinformatics_tools/tools/bcftools-1.18.sif bcftools view -R "$bed_file" "$cnvkit_vcf" -Oz -o "$subset_cnvkit_vcf"

intermediate_output_path="${output_dir}/${sample_id}.filtered.tsv"
python3 ./sv_insdeldup.py "$subset_manta_vcf" "$sample_id" "$subset_cnvkit_vcf" "$intermediate_output_path"

final_output_path="${output_dir}/${sample_id}_final.sv.tsv"
python3 ./validation_with_cnvkit.py "$intermediate_output_path" "$cnvkit_vcf" "$final_output_path"

# Remove intermediate files
rm $subset_manta_vcf $subset_cnvkit_vcf $intermediate_output_path
