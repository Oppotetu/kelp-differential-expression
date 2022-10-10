#!/bin/bash

#SBATCH -J move_align
#SBATCH --ntasks=30 # number CPUs to be used
#SBATCH --nodes=1
#SBATCH --mem=50G
#SBATCH --partition=smallmem
#SBATCH -o move_align

cd /mnt/SCRATCH/simonha/alignment

find . -name 'K*.bam' -exec mv -t ../corset/ {} \;

###
