#!/bin/bash

#SBATCH -J corset-reads_edit
#SBATCH --ntasks=150 # number CPUs to be used
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --partition=hugemem
#SBATCH -o corset-reads_slurm_edit

cd /mnt/SCRATCH/simonha/corset

####corset command###


###
# /mnt/users/simonha/software/corset-1.09-linux64/corset *.bam -r true-stop

for FILE in `ls *.bam` ; do
   /mnt/users/simonha/software/corset-1.09-linux64/corset -r true-stop $FILE &
done
wait
