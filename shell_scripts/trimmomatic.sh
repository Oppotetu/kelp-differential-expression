#!/bin/bash

## run with sbatch job_script

# Job name:
#SBATCH -J trimmomatic
# number CPUs to be used
#SBATCH -n 10
#SBATCH -N 1
#SBATCH --mem=10G
#SBATCH --partition=smallmem
#SBATCH -o trimmomatic_slurm

source activate qc

cd /path/to/raw/fq files
#
trimmomatic PE -threads 80 -phred33 S1_1.fq.gz S2_2.fq.gz \
                ../cleaned_data/S1_1.trimm.fq.gz ../cleaned_data/S1_1.singles.fq.gz \
                ../cleaned_data/S1_2.trimm.fq.gz ../cleaned_data/S1_2.singles.fq.gz \
                SLIDINGWINDOW:4:25 MINLEN:35



rm *singles.fq.gz   #deleting  singles files

##saves the output files in the pwd


