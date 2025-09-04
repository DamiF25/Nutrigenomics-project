#!/bin/bash
#SBATCH --job-name=FastQC_Trim
#SBATCH --mail-user=dfamakinde01@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem-per-cpu=10G
#SBATCH --output=../Logs/FastQC_Trim_Just_30.log
#SBATCH --partition=bio-compute
#SBATCH --chdir=/mnt/scratch2/users/40319381/Nutrigenomics/Read_Files

module load apps/fastqc/0.11.8/noarch
module load apps/trimgalore/0.4.4/noarch

# Provide paths to both forward and reverse read files for paired-end trimming
trim_galore --fastqc -q 30 -length 50 --paired PN0551_0030_S30_L001_R1_001.fastq.gz PN0551_0030_S30_L001_R2_001.fastq.gz
