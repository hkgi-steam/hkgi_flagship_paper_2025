#!/usr/bin/bash
# Extract variants resided on the targeted genes, i.e. CACNA1S, RYR1 and ABCG2
export sample_info="/home/_shared/jscliu/project/2025/Flagship/reference/sample_info_annot.2024-11-12.csv"
export data_dir="/mnt/data/patient_data/volume1/genomic_data/latest"
export projectDir="/home/_shared/jscliu/project/2025/Flagship/analysis/secondary/PGx/01d.genotype.manual"
export gene_regions="${projectDir}/script/gene_regions.tsv"
export out_dir="${projectDir}/results"
export bcftools_sif="/mnt/bioinformatics_tools/tools/bcftools-1.18.sif"
mkdir -p $out_dir

extract_var() {
    entry=$1
    sre_id=`echo $entry | cut -d"," -f1`
    sre_sample_id=`echo $entry | cut -d"," -f2`
    smallvar_vcf="${data_dir}/${sre_id}/${sre_sample_id}/${sre_sample_id}.smallvar.vcf.gz"
    genes_tsv="${out_dir}/${sre_id}.genes.tsv"
    [ ! -f $smallvar_vcf ] && exit    # Exit if the smallvar VCF is not available

    # Extract variants in the genes and convert into TSV
    echo -e "CHROM\tPOS\tREF\tALT\tGT\tAD" > $genes_tsv    # header
    apptainer run $bcftools_sif bcftools query \
        -R $gene_regions \
        -f "%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t[%AD]" \
        $smallvar_vcf >> $genes_tsv
    [ $? -eq 0 ] && echo $sre_id
}

export -f extract_var
grep "^SRE" $sample_info | parallel -j20 extract_var
