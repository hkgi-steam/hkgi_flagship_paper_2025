#!/usr/bin/bash
set -eu
export inputDir="/home/_shared/project/allsample/vep"
export toolsDir="/mnt/bioinformatics_tools/tools"
export dominant_gene_panel="/home/_shared/jscliu/project/2025/Flagship/reference/acmg_sf3_2.ext200bp.bed"
export outDir="/home/_shared/jscliu/project/2025/Flagship/analysis/secondary/dominant/01.known_variants/results"

vep_extract_af() {
    in_vcf=$1
    
    # Define variables from input
    sample_id=`basename $in_vcf | cut -d"." -f1`
    out_vcf="${outDir}/${sample_id}.split.rare.vep.additionalFinding.vcf.gz"
    out_tsv="${outDir}/${sample_id}.split.rare.vep.additionalFinding.tsv"

    # Skip if the vep VCF is already processed
    if [[ -f $out_tsv ]]; then 
        echo "Found and cached $sample_id"
        return 0
    fi

    # Extract variants of additional finding genes
    apptainer exec ${toolsDir}/bcftools-1.18.sif bcftools view --threads 4 $1 \
        -R $dominant_gene_panel \
        -Oz -o $out_vcf
    apptainer exec ${toolsDir}/bcftools-1.18.sif bcftools index -t $out_vcf
    
    # Convert VCF to TSV
    echo -e "CHROM\tPOS\tREF\tALT\tGT\tAD\tGene\tHGVSc\tHGVSp\tConsequence\tREVEL\tSpliceAI\tclinvar.vcf.gz\tclinvar.vcf.gz_CLNSIG\tclinvar.vcf.gz_CLNREVSTAT" > $out_tsv
    apptainer exec ${toolsDir}/bcftools-1.18.sif bcftools +split-vep $in_vcf \
        -f "%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t[%AD]\t%SYMBOL\t%HGVSc\t%HGVSp\t%Consequence\t%REVEL\t%spliceai_snv_hg38_gt02.vcf.gz_SpliceAI\t%clinvar.vcf.gz\t%clinvar.vcf.gz_CLNSIG\t%clinvar.vcf.gz_CLNREVSTAT" >> $out_tsv
}

export -f vep_extract_af
ls ${inputDir}/*.split.rare.vep.vcf.gz | parallel -j12 vep_extract_af
