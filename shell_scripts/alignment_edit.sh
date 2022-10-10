#!/bin/bash

## run with sbatch job_script
# Job name:
#SBATCH -J RSEM_edit
#SBATCH --ntasks=130 # number CPUs to be used
#SBATCH --nodes=1
#SBATCH --mem=130G
#SBATCH --partition=hugemem
#SBATCH -o alignment_edit_slurm
#SBATCH --error alignment_edit_error

cd /mnt/SCRATCH/simonha/alignment

#quantification
for infile in K1_1.clean.fq.gz K2_1.clean.fq.gz K3_1.clean.fq.gz K4_1.clean.fq.gz K5_1.clean.fq.gz K6_1.clean.fq.gz K7_1.clean.fq.gz K8_1.clean.fq.gz K10_1.clean.fq.gz K11_1.clean.fq.gz K12_1.clean.fq.gz K13_1.clean.fq.gz K14_1.clean.fq.gz K15_1.clean.fq.gz K16_1.clean.fq.gz K17_1.clean.fq.gz K18_1.clean.fq.gz K19_1.clean.fq.gz K20_1.clean.fq.gz K21_1.clean.fq.gz K22_1.clean.fq.gz K23_1.clean.fq.gz K24_1.clean.fq.gz K25_1.clean.fq.gz K26_1.clean.fq.gz K27_1.clean.fq.gz K28_1.clean.fq.gz K29_1.clean.fq.gz K30_1.clean.fq.gz K31_1.clean.fq.gz K32_1.clean.fq.gz K33_1.clean.fq.gz K34_1.clean.fq.gz K35_1.clean.fq.gz K36_1.clean.fq.gz K37_1.clean.fq.gz K38_1.clean.fq.gz K39_1.clean.fq.gz K40_1.clean.fq.gz K41_1.clean.fq.gz K42_1.clean.fq.gz K43_1.clean.fq.gz K44_1.clean.fq.gz;
do
base=$(basename ${infile} _1.clean.fq.gz)
align_and_estimate_abundance.pl --seqType fq --SS_lib_type FR \
                                --transcripts /mnt/SCRATCH/simonha/alignment/unseqs.fasta \
                                --left <(zcat ${infile}) --right <(zcat ${base}_2.clean.fq.gz) --est_method RSEM --aln_method bowtie2 \
                                --trinity_mode --prep_reference --output_dir /mnt/SCRATCH/simonha/alignment/${base} --thread_count 80 2>&1 > ${base}_opt
done

#S1_1.trimm.fq.gz  S2_1.trimm.fq.gz  S49_1.trimm.fq.gz  S38_1.trimm.fq.gz S43_1.trimm.fq.gz  S44_1.trimm.fq.gz  S7_1.trimm.fq.gz  S8_1.trimm.fq.gz
