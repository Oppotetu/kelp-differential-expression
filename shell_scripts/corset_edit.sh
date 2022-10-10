#!/bin/bash

## run with sbatch job_script

# Job name:
#SBATCH -J corset_edit
#SBATCH --ntasks=50 # number CPUs to be used
#SBATCH --nodes=1
#SBATCH --mem=50G
#SBATCH --partition=hugemem
#SBATCH -o corset_slurm_edit

cd /mnt/SCRATCH/simonha/corset

####corset command###

####

corset K1.bam.corset-reads K2.bam.corset-reads K3.bam.corset-reads K4.bam.corset-reads K5.bam.corset-reads K6.bam.corset-reads \
K7.bam.corset-reads K8.bam.corset-reads K10.bam.corset-reads K11.bam.corset-reads K12.bam.corset-reads K13.bam.corset-reads \
K14.bam.corset-reads K15.bam.corset-reads K16.bam.corset-reads K17.bam.corset-reads K18.bam.corset-reads K19.bam.corset-reads \
K20.bam.corset-reads K21.bam.corset-reads K22.bam.corset-reads K23.bam.corset-reads K24.bam.corset-reads K25.bam.corset-reads \
K26.bam.corset-reads K27.bam.corset-reads K28.bam.corset-reads K29.bam.corset-reads K30.bam.corset-reads K31.bam.corset-reads \
K32.bam.corset-reads K33.bam.corset-reads K34.bam.corset-reads K35.bam.corset-reads K36.bam.corset-reads K37.bam.corset-reads \
K38.bam.corset-reads K39.bam.corset-reads K40.bam.corset-reads K41.bam.corset-reads K42.bam.corset-reads K43.bam.corset-reads \
K44.bam.corset-reads \
 -n K1,K2,K3,K4,K5,K6,K7,K8,K10,K11,K12,K13,K14,K15,K16,K17,K18,K19,K20,K21,K22,K23,K24,K25,K26,K27,K28,K29,K30,K31,K32,K33,K34,K35,K36,K37,K38,K39,K40,K41,K42,K43,K44 \
 -g 1,1,1,1,2,2,2,2,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6,7,7,7,7,8,8,8,8,9,9,9,9,10,10,10,10,11,11,11,11 \
 -d 0.3 -m 10 -i corset -f true \

