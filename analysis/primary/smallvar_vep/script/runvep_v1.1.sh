invcf=$1
outvcf=$2

apptainer run \
  --bind /home/:/home/ \
  --bind /mnt/:/mnt/ \
  /mnt/bioinformatics_tools/tools/vep-110.0.sif vep \
  --offline \
  --cache \
  --dir_cache /mnt/bioinformatics_tools/data/vep/ \
  --dir_plugins /mnt/bioinformatics_tools/data/other/vep_plugins/vep_plugins_perl \
  --species homo_sapiens \
  --refseq \
  --assembly GRCh38 \
  --fasta /mnt/bioinformatics_tools/data/Homo_sapiens_assembly38.fasta \
  --pick \
  --pick_order rank,mane_select \
  --use_given_ref \
  --hgvs \
  --symbol \
  --force_overwrite \
  --check_existing \
  --clin_sig_allele 1 \
  -i ${invcf} \
  --format vcf \
  -o ${outvcf} \
  --stats_text \
  --vcf \
  --plugin REVEL,/mnt/bioinformatics_tools/data/other/vep_plugins/new_tabbed_revel_grch38.tsv.gz,1 \
  --custom file=/home/_shared/database/gnomad/gnomad.genomes.v3.1.2.sites.allchr_eas_unsort.vcf.gz,short_name=gnomad,format=vcf,type=exact,fields=AC%AN%AF%AC_eas%AN_eas%AF_eas \
  --custom file=/home/_shared/database/vep/spliceai_snv_hg38_gt02.vcf.gz,shortname=spliceai,format=vcf,type=exact,fields=SpliceAI \
  --custom file=/home/_shared/database/vep/clinvar/clinvar.vcf.gz,shortname=ClinVar,format=vcf,type=exact,fields=CLNSIG%CLNREVSTAT%CLNDN

