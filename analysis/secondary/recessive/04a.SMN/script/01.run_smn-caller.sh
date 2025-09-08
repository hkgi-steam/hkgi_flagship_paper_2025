#!/bin/bash
manifest="full_cohort_sr_path.txt"
reference="/home/_shared/database/reference/Homo_sapiens_assembly38.fasta"
out_dir="/home/_shared/jscliu/project/2025/Flagship/analysis/secondary/recessive/04a.SMN/summary"

smn_caller_script="/SMNCopyNumberCaller-master/smn_caller.py"
prefix="smn_output"
threads=64
genome="38"

apptainer exec /home/_shared/tools/smn-caller.sif \
    python3 $smn_caller_script --manifest $manifest \
                               --genome $genome \
                               --prefix $prefix \
                               --outDir $out_dir \
                               --threads $threads \
                               --reference $reference
