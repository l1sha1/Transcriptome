#!/bin/bash

#SBATCH --mem 150GB
#SBATCH -p  amd_256M
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 40:00:00


files=`cat files.txt` 
for NAME in $files
do
echo "$files"
gunzip ./${NAME}_R1_001.fastq.gz
gunzip ./${NAME}_R2_001.fastq.gz
done

