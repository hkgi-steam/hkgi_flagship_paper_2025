#!/usr/bin/bash
projectDir="/home/_shared/jscliu/project/2025/Flagship/analysis/secondary/recessive/04b.HBA"
HBA_caller="${projectDir}/script/HBA_caller"
sample_bam_list="${projectDir}/script/samples_cram_path.txt"
tmpDir="${HBA_caller}/tmp"
outDir="${projectDir}/summary"

# Run the HBA caller
[[ -d $tmpDir ]] && rm -rf $tmpDir && mkdir $tmpDir
cd $HBA_caller
python3 main.py \
    --bam_list $sample_bam_list \
    --temp_dir $tmpDir \
    --output_dir $outDir \
    --threads 8 --chunk_size 100 
