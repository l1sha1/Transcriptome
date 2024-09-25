#!/bin/bash

#SBATCH --mem 250GB
#SBATCH -p  amd_1Tb
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 120:00:00


files=`cat files.txt` 

for NAME in $files

do
STAR --runThreadN 8 --runMode alignReads \
--genomeDir index \
--readFilesIn ./libs/${NAME}_R1_001.fastq ./libs/${NAME}_R2_001.fastq \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix ./$NAME \

done

