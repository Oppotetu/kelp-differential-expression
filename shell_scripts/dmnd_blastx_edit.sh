#!/bin/bash

# Job name:
#SBATCH -J blastx_dmnd
# number CPUs to be used
#SBATCH -n 25
#SBATCH -N 1
#SBATCH --mem=25G
#SBATCH --partition=smallmem
#SBATCH -o slurm_blastx_dmnd

cd /mnt/SCRATCH/simonha/diamond

diamond blastx --db ./phaeophyceae_db.dmnd -q ./light1_extracted.fasta \
--out ./output/.tsv --outfmt 6 qseqid sseqid bitscore evalue pident salltitles \
--ultra-sensitive --threads 88 --max-target-seqs 5

