#!/bin/bash

# Job name:
#SBATCH -J eggnog
#SBATCH -n 40 # number CPUs to be used
#SBATCH -N 1
#SBATCH --mem=30G
#SBATCH --partition=smallmem
#SBATCH -o slurm_eggnog

conda activate blast_apps

cd /net/cn-1/mnt/SCRATCH/IPV-RNAseq22/files

#create_dbs.py -m diamond --dbname streptophyta --taxa  --taxids 35493 --data_dir /mnt/SCRATCH2/IPV-RNAseq22/DBs/eggnog_db -y 

emapper.py -m diamond --itype CDS -i trinity_decontaminated1.Trinity.fasta --data_dir /mnt/SCRATCH2/IPV-RNAseq22/DBs/eggnog_db \
 -o Tim --output_dir ./ --excel --cpu 120 --dmnd_db /mnt/SCRATCH2/IPV-RNAseq22/DBs/eggnog_db/streptophyta.dmnd --usemem --dbtype seqdb --qtype seq \
 --tax_scope 35493 --target_taxa 4565,4513,15368,200361,4567,4572 --tax_scope_mode narrowest --override --report_orthologs --target_orthologs one2one --dbmem --sensmode ultra-sensitive



conda deactivate

#4565 wheat, 4513 barley, 15368 brachypodium, 200361 aegilops tauschi strangulata, 4567 triticum turgidum durum, 4572 triticum untaru