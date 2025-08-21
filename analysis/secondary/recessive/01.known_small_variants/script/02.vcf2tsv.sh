#!/usr/bin/bash
set -eu
export resultDir="/home/_shared/jscliu/project/2025/Flagship/analysis/secondary/recessive/01.known_variants/results"
export clinvar_vcf="/home/_shared/jscliu/project/2025/Flagship/analysis/secondary/recessive/01.known_variants/data/recessive.known_variants.hgnc.vf.gz"
export toolsDir="/mnt/bioinformatics_tools/tools"

vcf2tsv() {
    vcf=$1
    tsv=`echo $vcf | sed "s/\.vcf\.gz/\.tsv/g"`

    # Define columns for bcftools query
    if [[ $(basename $vcf) =~ ^clinvar* ]]; then 
        cols="%CHROM\t%POS\t%REF\t%ALT\t%ID\t%INFO/GENEINFO\t%INFO/CLNSIG\t%INFO/CLNREVSTAT"
    elif [[ $(basename $vcf) =~ ^SRE.* ]]; then 
        cols="%CHROM\t%POS\t%REF\t%ALT\t[%GT]"
    fi

    # bcftools norm then query
    if [[ -f $tsv ]]; then 
        echo "Cached: $tsv"
    else
        apptainer exec ${toolsDir}/bcftools-1.18.sif bcftools norm \
            -m-any $vcf -Ov | \
            apptainer exec ${toolsDir}/bcftools-1.18.sif bcftools query \
            -f $cols - > $tsv
        if [[ $? -eq 0 ]]; then echo "Completed: ${vcf}"; else echo "Failed: ${vcf}"; fi
    fi
}

export -f vcf2tsv
vcf2tsv $clinvar_vcf
find $resultDir -maxdepth 1 -type f -name "*.smallvar.recessive.vcf.gz" | parallel -j24 vcf2tsv
