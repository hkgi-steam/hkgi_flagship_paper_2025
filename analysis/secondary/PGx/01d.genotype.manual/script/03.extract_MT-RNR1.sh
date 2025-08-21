#!/usr/bin/bash
set -eu

# Extract variants resided on MT-RNR1
export sample_info="/home/_shared/jscliu/project/2025/Flagship/reference/sample_info_annot.2024-11-12.csv"
export data_dir="/mnt/data/patient_data/volume1/genomic_data/latest"
export projectDir="/home/_shared/jscliu/project/2025/Flagship/analysis/secondary/PGx/01d.genotype.manual"
export out_dir="${projectDir}/results"
export bcftools_sif="/mnt/bioinformatics_tools/tools/bcftools-1.18.sif"
mkdir -p $out_dir

extract_mt_var() {
    entry=$1
    sre_id=`echo $entry | cut -d"," -f1`
    sre_sample_id=`echo $entry | cut -d"," -f2`
    mito_vcf="${data_dir}/${sre_id}/${sre_sample_id}/${sre_sample_id}.mito.vcf"
    tmp_vcf="${out_dir}/${sre_id}.mito.vcf.gz"
    genes_tsv="${out_dir}/${sre_id}.MT-RNR1.tsv"
    [ ! -f $mito_vcf ] && exit    # Exit if the smallvar VCF is not accessible

    # Create a copy of VCF and bgzip
    apptainer run $bcftools_sif bcftools view \
        --write-index=tbi -Oz -o $tmp_vcf $mito_vcf

    # Extract variants in MT-RNR1 and convert into TSV
    echo -e "CHROM\tPOS\tREF\tALT\tGT\tAD\tAF" > $genes_tsv    # header
    apptainer run $bcftools_sif bcftools query \
        -r chrM:648-1601 \
        -f "%CHROM\t%POS\t%REF\t%ALT\t[ %GT]\t[%AD]\t[%AF]" \
        $tmp_vcf >> $genes_tsv
    
    [ $? -eq 0 ] && echo $sre_id && rm $tmp_vcf && rm ${tmp_vcf}.tbi
}

export -f extract_mt_var
grep "^SRE" $sample_info | parallel -j20 extract_mt_var
