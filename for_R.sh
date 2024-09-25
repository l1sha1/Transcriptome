#!/bin/bash

#SBATCH --mem 100GB
#SBATCH -p  amd_256M
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 100:00:00

export PATH="/home/lishaiea/miniconda3/bin:$PATH"
source activate rna-seq

cd /beegfs/scratch/ws/ws1/lishaiea-work/

./Deseq2.R
