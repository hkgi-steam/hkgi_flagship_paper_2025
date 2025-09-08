#!/usr/bin/bash
set -eu
export clinvar_vcf="/home/_shared/jscliu/project/2025/Flagship/analysis/secondary/recessive/01.known_small_variants/data/recessive.known_variants.hgnc.vf.gz"
export toolsDir="/mnt/bioinformatics_tools/tools"
export cohort_list="/home/_shared/jscliu/project/2025/Flagship/reference/sample_info_annot.2024-11-12.csv"
export dataDir="/mnt/data/patient_data/volume1/genomic_data/latest"
export resultDir="/home/_shared/jscliu/project/2025/Flagship/analysis/secondary/recessive/01.known_small_variants/results"

extract_var() {
    cohort_individual=$1
    sre_participant_id=`echo $cohort_individual | cut -d"," -f1`
    sre_sample_id=`echo $cohort_individual | cut -d"," -f2`
    
    smallvar_vcf="${dataDir}/${sre_participant_id}/${sre_sample_id}/${sre_sample_id}.smallvar.vcf.gz"
    out_vcf="${resultDir}/${sre_sample_id}.smallvar.recessive.vcf.gz"
    
     # Skip when the index file is already there
    if [[ -f ${out_vcf}.tbi ]]; then
        echo "Found and cached: $sre_sample_id"
    else
        # bcftools view to subset
        apptainer exec ${toolsDir}/bcftools-1.18.sif bcftools view --threads 5 $smallvar_vcf \
            -R $clinvar_vcf \
            -Oz -o $out_vcf
        
        # bcftools index
        apptainer exec ${toolsDir}/bcftools-1.18.sif bcftools index -t $out_vcf  

        if [ $? -eq 0 ]; then 
            echo "Done: $sre_sample_id"
        fi
    fi
}

export -f extract_var
tail -n +2 $cohort_list | parallel -j12 extract_var
