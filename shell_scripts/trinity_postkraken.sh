#!/bin/bash

## run with sbatch job_script

# Job name:
#SBATCH -J postkraken_trinity
#SBATCH -n 130 # number CPUs to be used
#SBATCH -N 1
#SBATCH --mem=110G
#SBATCH --partition=hugemem
#SBATCH -o trinity_slurm_postkraken3
#SBATCH --error trinity_pdb_error
#SBATCH --mail-user=simon.oscar.andersen.haughom@gmail.com
#SBATCH --mail-type=ALL

cd /mnt/SCRATCH/simonha/postkraken_trinity_assembly

# trinity denovo assembly

Trinity --seqType fq --max_memory 110G --CPU 110 \
           --left uK1_1.fq uK2_1.fq uK10_1.fq uK11_1.fq uK13_1.fq uK14_1.fq uK17_1.fq uK18_1.fq \
           --right uK1_2.fq uK2_2.fq uK10_2.fq uK11_2.fq uK13_2.fq uK14_2.fq uK17_2.fq uK18_2.fq \
           --min_contig_length 200 --SS_lib_type FR --full_cleanup \
           --min_kmer_cov 2 --output /mnt/SCRATCH/simonha/postkraken_trinity_assembly \


