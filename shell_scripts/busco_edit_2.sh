#!/bin/bash

## run with sbatch job_script

#SBATCH -J busco_edit_2
#SBATCH --ntasks=30 # number CPUs to be used
#SBATCH --nodes=1
#SBATCH --mem=30G
#SBATCH --partition=smallmem
#SBATCH -o busco_slurm_2


cd /mnt/SCRATCH/simonha/

busco -m transcriptome -i ./safekelp_trinity_assembly.Trinity.fasta -o busco_4_output -l stramenopiles_odb10 -c 30 

# -f

#busco command
#-m denotes method to choose
#-i denotes the input file / specify the path to inputfile if it is not in current working directory
#-o denotes the output file / it gives a folder with the specified name  in current pwd
#-l denotes the lineage to be used /viridiplanate_odb10 indicates the plants lineage in BUSCO database
#-f denotes to run the analysis by force /useful
#-c cpus


#output it will be in the  pwd. Do not provide a path to save output in a directory
