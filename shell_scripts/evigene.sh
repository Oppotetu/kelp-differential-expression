#!/bin/bash

## run with sbatch job_script

# Job name:
#SBATCH -J evigene
#SBATCH --ntasks=30 # number CPUs to be used
#SBATCH --nodes=1
#SBATCH --mem=60G
#SBATCH --partition=hugemem
#SBATCH -o evigene_slurm
#SBATCH --error evigene_slurm_error


cd /mnt/SCRATCH/simonha

####corset command###

####

 /mnt/users/simonha/software/evigene/scripts/prot/tr2aacds.pl -log \
 -cdna /mnt/SCRATCH/simonha/safekelp_trinity_assembly.Trinity.fasta -NCPU 30 -MAXMEM 60000 \
