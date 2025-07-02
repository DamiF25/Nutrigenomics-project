#!/bin/bash
#SBATCH --job-name=Stringtie
#SBATCH --mail-user=email
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem-per-cpu=10G
#SBATCH --array=1-17
#SBATCH --output=./Logs/SAM_BAM_%j.log
#SBATCH --partition=bio-compute
#SBATCH --chdir=/path/to/file

module load apps/samtools/1.9/gcc-4.8.5

# Specify the path to the config file
config=/path/to/myarray.txt

# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

echo "This is array task ${SLURM_ARRAY_TASK_ID}, the sample name is ${sample}"

samtools sort -@ 4 -o ${sample}_Bgla_aligned.bam ${sample}_aligned.sam
