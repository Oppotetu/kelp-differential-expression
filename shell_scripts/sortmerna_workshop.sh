#!/bin/bash

## run with sbatch job_script

# Job name:
#SBATCH -J sortmerna
#SBATCH -n 80 # number CPUs to be used
#SBATCH -N 1
#SBATCH --mem=100G
#SBATCH --partition=hugemem
#SBATCH -o sortmerna_slurm

source activate decontamination

cd /path/to working directory

#
sortmerna --ref /mnt/SCRATCH/IPV-RNAseq22/DBs/sortmerna/rfam-5.8s-database-id98.fasta \
  --ref /mnt/SCRATCH/IPV-RNAseq22/DBs/sortmerna/rfam-5s-database-id98.fasta \
  --ref /mnt/SCRATCH/IPV-RNAseq22/DBs/sortmerna/silva-arc-16s-id95.fasta \
  --ref /mnt/SCRATCH/IPV-RNAseq22/DBs/sortmerna/silva-arc-23s-id98.fasta \
  --ref /mnt/SCRATCH/IPV-RNAseq22/DBs/sortmerna/silva-bac-16s-id90.fasta \
  --ref /mnt/SCRATCH/IPV-RNAseq22/DBs/sortmerna/silva-bac-23s-id98.fasta \
  --ref /mnt/SCRATCH/IPV-RNAseq22/DBs/sortmerna/silva-euk-18s-id95.fasta \
  --ref /mnt/SCRATCH/IPV-RNAseq22/DBs/sortmerna/silva-euk-28s-id98.fasta \
  --workdir /path/to/prefered/working directory/ \
  --idx-dir /mnt/SCRATCH/IPV-RNAseq22/DBs/sortmerna/index \
  --reads /path/to/input/filename.fasta \
  --fastx --blast 1 --num_alignments 1 --threads 80 --v
