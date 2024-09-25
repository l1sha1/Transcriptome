#!/bin/bash

#SBATCH --mem 100GB
#SBATCH -p  amd_256M
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 5:00:00


STAR --runThreadN 4 --runMode genomeGenerate --genomeDir index \
--genomeFastaFiles ./GCF_017639785.1_BCM_Maur_2.0_genomic.fna --sjdbGTFfile ./GCF_017639785.1_BCM_Maur_2.0_genomic.gff \
--sjdbOverhang 99

