#!/bin/bash
i_dir="/mnt/data/patient_data/volume1/genomic_data/latest"
o_dir="/home/_shared/jscliu/project/2025/Flagship/analysis/primary/01.smalllvar_vep/results"

prefix=$1

mkdir -p "$o_dir"
for subdir in "$i_dir"/*; do
  if [[ -d "$subdir" ]]; then
    sample_id=$(basename "$subdir")
    if [[ $sample_id == $prefix* ]]; then
      echo "$sample_id"
      i_vcf="$subdir/${sample_id}B01/${sample_id}B01.smallvar.vcf.gz"
      split_vcf="$o_dir/$sample_id.split.vcf.gz"
      split_rare_vcf="$o_dir/$sample_id.split.rare.vcf.gz"
      split_rare_vep_vcf="$o_dir/$sample_id.split.rare.vep.vcf"
      apptainer run /home/scyy/Data/bioinformatics_tools/tools/bcftools-1.18.sif bcftools norm -m -any -f ~/Data/bioinformatics_tools/data/Homo_sapiens_assembly38.fasta "$i_vcf" -Oz -o "$split_vcf"
      tabix -p vcf "$split_vcf"
      apptainer run /home/scyy/Data/bioinformatics_tools/tools/bcftools-1.18.sif bcftools isec -C -w1 "$split_vcf" /home/_shared/database/gnomad/gnomad.genomes.v3.1.2.sites.allchr_afgt0.005.vcf.gz -Oz -o "$split_rare_vcf"
      tabix -p vcf "$split_rare_vcf" 
      bash ./runvep_v1.1.sh "$split_rare_vcf" "$split_rare_vep_vcf"
      bgzip "$split_rare_vep_vcf"
      tabix -p vcf "$split_rare_vep_vcf.gz"
      rm "$split_vcf"
      rm "$split_vcf.tbi"
      rm "$split_rare_vcf"
      rm "$split_rare_vcf.tbi"
    fi
  fi
done
