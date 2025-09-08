#!/usr/bin/bash
set -eu

aldy_sif="/home/_shared/tools/aldy-v4.6.sif"
data_dir="/mnt/data/patient_data/volume1/genomic_data/latest"
ref_genome="/mnt/bioinformatics_tools/data/Homo_sapiens_assembly38.fasta"
outDir="/home/_shared/jscliu/project/2025/Flagship/analysis/secondary/PGx/01a.genotype.aldy/results"
sre_list=$1

cat $sre_list | parallel -j12 --line-buffer '
sre_participant_id={1}
data_dir={2}
aldy_sif={3}
ref_genome={4}
outDir={5}

cram=$(ls ${data_dir}/${sre_participant_id}/*/*.cram)
sample_id=`basename $cram | cut -d"." -f1`

if [ -f ${outDir}/${sample_id}.aldy ]; then
    echo "[LOG] Found and cached: ${sample_id}"
else
    apptainer run $aldy_sif aldy genotype \
        --profile wgs \
        --reference $ref_genome \
        --output ${outDir}/tmp_${sample_id}.aldy \
        $cram
    [[ $? -eq 0 ]] && mv ${outDir}/tmp_${sample_id}.aldy ${outDir}/${sample_id}.aldy
fi
' :::: - ::: $data_dir ::: $aldy_sif ::: $ref_genome ::: $outDir
