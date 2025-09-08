#!/usr/bin/bash
export cohort_list="/home/_shared/jscliu/project/2025/Flagship/reference/sample_info_annot.2024-11-12.csv"
export data_dir="/mnt/data/patient_data/volume1/genomic_data/latest"
export projectDir="/home/_shared/jscliu/project/2025/Flagship/analysis/secondary/PGx/01d.genotype.manual"
export out_tsv="${projectDir}/results/IFNL3_4.consolidate.tsv"
export bcftools_sif="/mnt/bioinformatics_tools/tools/bcftools-1.18.sif" 

extract_IFNL() {
    entry=$1
    participant_id=`echo $entry | cut -d"," -f1`
    sample_id=`echo $entry | cut -d"," -f2`
    smallvar_vcf="${data_dir}/${participant_id}/${sample_id}/${sample_id}.smallvar.vcf.gz"

    apptainer run $bcftools_sif bcftools view \
        -r chr19:39252525,chr19:39248513,chr19:39248515 $smallvar_vcf | \
        apptainer run $bcftools_sif bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]' | \
        awk -v participant_id="$participant_id" '{print participant_id"\t"$0}'
}

export -f extract_IFNL
echo -e "participant_id\tchr\tpos\tref\talt\tGT" > $out_tsv
cat $cohort_list | parallel extract_IFNL >> $out_tsv
