#!/usr/bin/bash
export inputDir="/home/_shared/project/allsample/vep"
export gene_panel="/home/_shared/jscliu/project/2025/Flagship/reference/hkgi_recessive.ext200bp.bed"
export outDir="/home/_shared/jscliu/project/2025/Flagship/analysis/secondary/recessive/02.novel_small_variants/results"
export toolsDir="/mnt/bioinformatics_tools/tools"
mkdir -p $outDir

extract_novel_recessive() {
    # Define variables from input
    sample_id=`basename $1 | cut -d"." -f1`
    out_tsv="${outDir}/${sample_id}.split.rare.vep.recessive.tsv"

    # Skip if the vep VCF is already processed
    if [[ -f $out_tsv ]]; then 
        echo "Found and cached $sample_id"
        return 0
    fi

    # Extract variants of recessive genes as TSV
    echo -e "CHROM\tPOS\tREF\tALT\tGT\tAD\tGene\tHGVSc\tHGVSp\tConsequence\tREVEL\tSpliceAI\tclinvar.vcf.gz\tclinvar.vcf.gz_CLNSIG\tclinvar.vcf.gz_CLNREVSTAT" > $out_tsv
    apptainer exec ${toolsDir}/bcftools-1.18.sif bcftools view --threads 8 $1 -R $gene_panel | \
        apptainer exec ${toolsDir}/bcftools-1.18.sif bcftools +split-vep - \
        -f "%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t[%AD]\t%SYMBOL\t%HGVSc\t%HGVSp\t%Consequence\t%REVEL\t%spliceai_snv_hg38_gt02.vcf.gz_SpliceAI\t%clinvar.vcf.gz\t%clinvar.vcf.gz_CLNSIG\t%clinvar.vcf.gz_CLNREVSTAT" | \
        awk -F'\t' '$13 == "."' >> $out_tsv
}

export -f extract_novel_recessive
ls ${inputDir}/*.split.rare.vep.vcf.gz | parallel -j5 extract_novel_recessive
