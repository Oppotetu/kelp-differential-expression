#!/bin/bash

## run with sbatch job_script

# Job name:
#SBATCH -J Trinity
#
#SBATCH -n 110 # number CPUs to be used
#SBATCH -N 1
#SBATCH --mem=110G
#SBATCH --partition=hugemem
#SBATCH -o trinity_slurm_2
#SBATCH --mail-user=simon.oscar.andersen.haughom@gmail.com
#SBATCH --mail-type=ALL

module load Anaconda3/5.3.0

source activate Trinity

cd /mnt/SCRATCH/simonha/link_trinity2

# trinity denovo assembly

Trinity --seqType fq --max_memory 110G --CPU 110 \
           --left forward_1.fq.gz --right reverse_2.fq.gz \
           --min_contig_length 200 --SS_lib_type FR --full_cleanup \
           --min_kmer_cov 2 --output /mnt/SCRATCH/simonha/safekelp_trinity_assembly 

source deactivate