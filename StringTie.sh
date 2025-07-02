#!/bin/bash
#SBATCH --job-name=Stringtie
#SBATCH --mail-user=email
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem-per-cpu=10G
#SBATCH --array=1-17
#SBATCH --output=./Logs/Stringtie_%j.log
#SBATCH --partition=bio-compute
#SBATCH --chdir=/path/to/file

module load stringtie/1.3.6
module load apps/python/2.7.8/gcc-4.8.5

# Specify the path to the config file
config=/path/to/myarray.txt

# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

stringtie -e -B -G ./Genome/Bgla_Geno_Annos_No_Gene_Processed.gtf -o ./${sample}.gtf -p 4 -A ${sample}_Counts ${sample}_Bgla_aligned.bam
