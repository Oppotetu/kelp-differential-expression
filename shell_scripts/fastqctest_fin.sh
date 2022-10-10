#!/bin/bash

## run with sbatch job_script

# Job name:
#SBATCH -J fastqc
#
# number CPUs to be used
#SBATCH -n 20
#SBATCH -N 1
#SBATCH --mem=15G  
#SBATCH -o slurm_fastqc.out

module load Anaconda3/5.3.0

source activate Quality

cd /mnt/users/simonha/RNA_SEQ/link

# fastqc

fastqc *.fq.gz

source deactivate