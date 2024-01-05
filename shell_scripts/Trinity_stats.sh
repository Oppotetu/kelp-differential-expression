#!/bin/bash

## run with sbatch job_script

# Job name:
#SBATCH -J Trinity
#
#SBATCH -n 2 # number CPUs to be used
#SBATCH -N 1
#SBATCH -o trinity_stats



cd /net/fs-2/scale/OrionStore/Home/simonha/Trinity_stats

# trinity stats

TrinityStats.pl ./safekelp_trinity_assembly.Trinity.fasta >& stats.log

