#!/usr/bin/bash
export sample_info_csv="/home/_shared/jscliu/project/2025/Flagship/reference/sample_info_annot.2024-11-12.csv"
export data_dir="/mnt/data/patient_data/volume1/genomic_data"


## Sort and gzip the VCF from previous jupyter notebook
cd /home/_shared/jscliu/project/2025/Flagship/analysis/secondary/dominant/03.filtering/results
ls *.vcf | parallel -j 25 '
apptainer run /mnt/bioinformatics_tools/tools/bcftools-1.18.sif bcftools sort -Oz -o {}.gz {}
'
rm *.vcf

## Retrieve original INFO field for hard filtering
query_var() {
    # Define variables
    in_vcf=$1
    participant_id=`basename $in_vcf | cut -d"." -f1`
    sample_id=`grep -w $participant_id $sample_info_csv | cut -d"," -f2`
    data_vcf="${data_dir}/${participant_id}/${sample_id}/${sample_id}.smallvar.vcf.gz"
    out_vcf="${participant_id}.gatk_info.vcf.gz"
    
    # bcftools view to subset
    apptainer run /mnt/bioinformatics_tools/tools/bcftools-1.18.sif bcftools view \
        $data_vcf -R $in_vcf | \
        apptainer run /mnt/bioinformatics_tools/tools/bcftools-1.18.sif bcftools norm \
        -Oz -o $out_vcf -m - -

    if [[ $? -eq 0 ]]; then echo "[LOG] Done - ${participant_id}"; fi    # Logging
}

export -f query_var
ls *.query.vcf.gz | parallel -j20 query_var
