#!/bin/bash

## run with sbatch job_script

# Job name:
#SBATCH -J abundance_est_matrix_edit
#SBATCH --ntasks=75 # number CPUs to be used
#SBATCH --nodes=1
#SBATCH --mem=250G
#SBATCH --partition=hugemem
#SBATCH -o abundance_matrix_slurm_edit


cd /mnt/SCRATCH/simonha/alignment

# use find . -type f -name "*.isoforms.results" to find paths to all files  

abundance_estimates_to_matrix.pl --est_method RSEM --out_prefix FT --cross_sample_norm TMM --name_sample_by_basedir \
--gene_trans_map /mnt/SCRATCH/simonha/safekelp_trinity_assembly.Trinity.fasta.gene_trans_map \
./K1/RSEM.isoforms.results \
./K2/RSEM.isoforms.results \
./K3/RSEM.isoforms.results \
./K4/RSEM.isoforms.results \
./K5/RSEM.isoforms.results \
./K6/RSEM.isoforms.results \
./K7/RSEM.isoforms.results \
./K8/RSEM.isoforms.results \
./K10/RSEM.isoforms.results \
./K11/RSEM.isoforms.results \
./K12/RSEM.isoforms.results \
./K13/RSEM.isoforms.results \
./K14/RSEM.isoforms.results \
./K15/RSEM.isoforms.results \
./K16/RSEM.isoforms.results \
./K17/RSEM.isoforms.results \
./K18/RSEM.isoforms.results \
./K19/RSEM.isoforms.results \
./K20/RSEM.isoforms.results \
./K21/RSEM.isoforms.results \
./K22/RSEM.isoforms.results \
./K23/RSEM.isoforms.results \
./K24/RSEM.isoforms.results \
./K25/RSEM.isoforms.results \
./K26/RSEM.isoforms.results \
./K27/RSEM.isoforms.results \
./K28/RSEM.isoforms.results \
./K29/RSEM.isoforms.results \
./K30/RSEM.isoforms.results \
./K31/RSEM.isoforms.results \
./K32/RSEM.isoforms.results \
./K33/RSEM.isoforms.results \
./K34/RSEM.isoforms.results \
./K35/RSEM.isoforms.results \
./K36/RSEM.isoforms.results \
./K37/RSEM.isoforms.results \
./K38/RSEM.isoforms.results \
./K39/RSEM.isoforms.results \
./K40/RSEM.isoforms.results \
./K41/RSEM.isoforms.results \
./K42/RSEM.isoforms.results \
./K43/RSEM.isoforms.results \
./K44/RSEM.isoforms.results \
