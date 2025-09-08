#!/usr/bin/bash
set -e

picard_sif=`ls /mnt/bioinformatics_tools/tools/* | grep picard | head -1`
hlahd_sif="/home/_shared/tools/hlahd-v1.7.0.sif"
REFERENCE="/mnt/bioinformatics_tools/data/Homo_sapiens_assembly38.fasta"
OUTDIR="/home/_shared/jscliu/project/2025/Flagship/analysis/secondary/PGx/01c.genotype.HLA-HD/results"
CRAM=$1
sampleID=`basename $CRAM | sed "s/\.cram$//"`
RUN_LOG="${OUTDIR}/${sampleID}/${sampleID}.run_hlahd.log"
SAMPLE_OUTDIR="${OUTDIR}/${sampleID}/${sampleID}.hlahd_out"
STATUS_LOG_PREFIX="${OUTDIR}/${sampleID}/${sampleID}.run_hlahd"

mkdir -p ${SAMPLE_OUTDIR}

##### 0. Skip sample if the OK log presents
if [[ -f ${STATUS_LOG_PREFIX}.ok ]]; then echo "${sampleID} has been processed" && exit 0; fi

script_start=`date +%s`

##### 1. Subset aligmments mapped to mgc region and convert to FASTQ
# 1.1 Extract reads mapped to the MHC region and their mates
samtools view -@ 8 -h \
    -T $REFERENCE $CRAM chr6:28,510,120-33,480,577 | \
    awk '{print $1}' | sort | uniq > \
    ${OUTDIR}/${sampleID}/${sampleID}.mhc.read_names.txt
samtools view -@ 8 -h -b \
    -N ${OUTDIR}/${sampleID}/${sampleID}.mhc.read_names.txt \
    -T $REFERENCE $CRAM > \
    ${OUTDIR}/${sampleID}/${sampleID}.mhc_and_mates.bam
# 1.2 Extract unmap reads
samtools view -@ 8 -f 4 \
    -T $REFERENCE $CRAM | \
    awk '{print $1}' | sort | uniq > \
    ${OUTDIR}/${sampleID}/${sampleID}.unmap.read_names.txt
samtools view -@ 8 -h \
    -N ${OUTDIR}/${sampleID}/${sampleID}.unmap.read_names.txt \
    -T $REFERENCE $CRAM > \
    ${OUTDIR}/${sampleID}/${sampleID}.unmapped_and_mates.bam
# 1.3 Merge bam files
[ -f ${OUTDIR}/${sampleID}/${sampleID}.merge.bam ] && rm ${OUTDIR}/${sampleID}/${sampleID}.merge.bam
samtools merge -@ 8 \
    -o ${OUTDIR}/${sampleID}/${sampleID}.merge.bam \
    ${OUTDIR}/${sampleID}/${sampleID}.mhc_and_mates.bam \
    ${OUTDIR}/${sampleID}/${sampleID}.unmapped_and_mates.bam > $RUN_LOG 2>&1    
# 1.4 Create FASTQ
[ -f ${OUTDIR}/${sampleID}/${sampleID}.hlatmp.1.fastq ] && rm ${OUTDIR}/${sampleID}/${sampleID}.*.fastq
apptainer exec $picard_sif picard SamToFastq -Xmx16G \
     I=${OUTDIR}/${sampleID}/${sampleID}.merge.bam \
     F=${OUTDIR}/${sampleID}/${sampleID}.hlatmp.1.fastq \
     F2=${OUTDIR}/${sampleID}/${sampleID}.hlatmp.2.fastq
# 1.5 Change fastq ID
cat ${OUTDIR}/${sampleID}/${sampleID}.hlatmp.1.fastq | \
    awk '{if(NR%4==1) $0=$0" 1"; print}' > \
    ${OUTDIR}/${sampleID}/${sampleID}.hla.1.fastq
cat ${OUTDIR}/${sampleID}/${sampleID}.hlatmp.2.fastq | \
    awk '{if(NR%4==1) $0=$0" 2"; print}' > \
    ${OUTDIR}/${sampleID}/${sampleID}.hla.2.fastq

CramToFastq_end=`date +%s`

##### 2. Run HLA-HD
apptainer exec $hlahd_sif hlahd.sh \
    -t 15 -f /hlahd.1.7.0/freq_data \
    ${OUTDIR}/${sampleID}/${sampleID}.hla.1.fastq \
    ${OUTDIR}/${sampleID}/${sampleID}.hla.2.fastq \
    /hlahd.1.7.0/HLA_gene.split.3.50.0.txt \
    /hlahd.1.7.0/dictionary \
    $sampleID \
    ${SAMPLE_OUTDIR} >> $RUN_LOG 2>&1

hlahd_end=`date +%s`
    
CramToFastq_timed_s=$(($CramToFastq_end-$script_start))
CramToFastq_timed=$(($CramToFastq_timed_s/60))
hlahd_timed_s=$(($hlahd_end-$CramToFastq_end))
hlahd_timed=$(($hlahd_timed_s))

##### 3. Clean-up temp_file and write to log
if [ $? -eq 0 ]; then 
    echo "[${sampleID}] took ${CramToFastq_timed} mins and ${hlahd_timed} secs to make FASTQ and run HLA-HD respectively" > ${STATUS_LOG_PREFIX}.ok
    rm ${OUTDIR}/${sampleID}/${sampleID}*.bam
    rm ${OUTDIR}/${sampleID}/${sampleID}*.fastq

    # Compress the SAM files in hlahd_out
    mapfile_dir="${SAMPLE_OUTDIR}/${sampleID}/mapfile"
    exon_dir="${SAMPLE_OUTDIR}/${sampleID}/exon"
    intron_dir="${SAMPLE_OUTDIR}/${sampleID}/intron"
    cwd=`pwd`
    [[ -d $mapfile_dir ]] && cd ${SAMPLE_OUTDIR}/${sampleID} && tar -czvf mapfile.tar.gz mapfile/ && rm -rf mapfile && cd $cwd
    [[ -d $exon_dir ]] && cd ${SAMPLE_OUTDIR}/${sampleID} && tar -czvf exon.tar.gz exon/ && rm -rf exon && cd $cwd
    [[ -d $intron_dir ]] && cd ${SAMPLE_OUTDIR}/${sampleID} && tar -czvf intron.tar.gz intron/ && rm -rf intron && cd $cwd
    
    # Remove read_names files
    [[ -f ${OUTDIR}/${sampleID}/${sampleID}.mhc.read_names.txt ]] && rm ${OUTDIR}/${sampleID}/${sampleID}.mhc.read_names.txt
    [[ -f ${OUTDIR}/${sampleID}/${sampleID}.unmap.read_names.txt ]] && rm ${OUTDIR}/${sampleID}/${sampleID}.unmap.read_names.txt
else
    touch ${STATUS_LOG_PREFIX}.err
fi

