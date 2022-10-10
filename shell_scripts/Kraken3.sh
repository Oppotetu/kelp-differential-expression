#!/bin/bash

## run with sbatch job_script

# Job name:
#SBATCH -J kraken_4
#SBATCH -n 110 # number CPUs to be used
#SBATCH -N 1
#SBATCH --mem=100G
#SBATCH --partition=hugemem
#SBATCH -o kraken_slurm_4

cd /mnt/SCRATCH/simonha/Kraken4

#to decontaminate a denovo assembly / single end reads

kraken2 -db /mnt/SCRATCH/simonha/Kraken2/db_deconta /mnt/SCRATCH/simonha/safekelp_trinity_assembly.Trinity.fasta \
--threads 70 --classified-out /mnt/SCRATCH/simonha/Kraken4/class/cseqs.fasta --unclassified-out /mnt/SCRATCH/simonha/Kraken4/unclass/unseqs.fasta \
--output . --report ./report_Kraken4

#paired end reads
# kraken2 --paired -db /mnt/SCRATCH/simonha/Kraken2/db_deconta /mnt/SCRATCH/simonha/safekelp_trinity_assembly.Trinity.fasta /mnt/SCRATCH/IPV-RNAseq22/cleaned_data/reads_2.fq.gz \
# --threads 130 --classified-out cseqs#.fq.gz --unclassified-out unseqs#.fq.gz \
# --output /mnt/SCRATCH/IPV-RNAseq22/kraken_opt/krkn_optput.tsv --report denovo_report --use-names  --confidence 0.05

##saves the output files in the pwd

