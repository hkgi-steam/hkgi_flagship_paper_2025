#!/usr/bin/bash
ls /mnt/data/patient_data/volume1/genomic_data/latest | parallel -j2 --line-buffer '
cram=$(ls /mnt/data/patient_data/volume1/genomic_data/latest/{}/{}*/{}*.cram)
bash run_hlahd.sh $cram
'