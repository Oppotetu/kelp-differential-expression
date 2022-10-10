#!/bin/bash

#SBATCH -J dmnd_make_db

# number CPUs to be used
#SBATCH -n 50
#SBATCH -N 1
#SBATCH --mem=60G
#SBATCH --partition=smallmem
#SBATCH -o slurm_dmnd_make

cd /mnt/SCRATCH/simonha/diamond

diamond makedb --in ./uniprot-taxonomy__Phaeophyceae.fasta -d ./phaeophyceae_db

#database file with extension databasename.dmnd with be output directory
