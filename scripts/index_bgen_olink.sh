#!/bin/bash
module load gcc/5.2.0
module load bgenix/1.0.1

# Create bgi index file for bgen files
cd /home/jm2294/GENETIC_DATA/INTERVAL/olink/bgen
CHR=$SLURM_ARRAY_TASK_ID
bgenix -g interval.imputed.olink.chr_${CHR}.bgen -index 