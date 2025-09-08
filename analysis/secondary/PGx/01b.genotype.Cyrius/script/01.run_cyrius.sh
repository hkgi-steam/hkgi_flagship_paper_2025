#!/usr/bin/bash
todays_date="2024-11-20"
indir_master="/mnt/data/patient_data/volume1/genomic_data"
projectDir="/home/_shared/jscliu/project/2025/Flagship/analysis/secondary/PGx/01b.genotype.Cyrius"
outdir="${projectDir}/results"
out_prefix="manifest.${todays_date}"
logdir="${projectDir}/log"
manifest="${projectDir}/script/manifest.${todays_date}.txt"
CYRIUS_SIF="/home/_shared/tools/cyrius-v1.1.1.sif"
GRCh38_fa="/mnt/bioinformatics_tools/data/Homo_sapiens_assembly38.fasta"
n_threads=36
mkdir -p $outdir

apptainer exec $CYRIUS_SIF python /usr/Cyrius/star_caller.py \
    --manifest $manifest \
    --genome 38 \
    --reference $GRCh38_fa \
    --outDir $outdir \
    --threads $n_threads \
    --prefix $out_prefix > ${logdir}/${out_prefix}.log 2>&1

