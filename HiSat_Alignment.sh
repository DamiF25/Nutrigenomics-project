#!/bin/bash 
#SBATCH --job-name=Alignment
#SBATCH --mail-user=email
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem-per-cpu=15G
#SBATCH --array=1-17
#SBATCH --output=./Logs/Alignment_%j.log
#SBATCH --partition=k2-gpu-v100
#SBATCH --chdir=/path/to/file

module load hisat2/2.1.0

# Specify the path to the config file
config=/path/to/myarray.txt

# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

# Paths to the forward and reverse read files
forward_read="/mnt/scratch2/users/40319381/New_Kraken/Read_Files/Read_Files/${sample}_R1_001_val_1.fq.gz"
reverse_read="/mnt/scratch2/users/40319381/New_Kraken/Read_Files/Read_Files/${sample}_R2_001_val_2.fq.gz"

echo "This is array task ${SLURM_ARRAY_TASK_ID}, the sample name is ${sample}"

hisat2 -p 4 -q --un-conc-gz ./${sample}_UNMAPPED -S ${sample}_aligned.sam -x ./Genomes/Bgla_Geno -1 $forward_read -2 $reverse_read
mv ${sample}_UNMAPPED.1 ${sample}_UNMAPPED.1.fq.gz
gunzip ${sample}_UNMAPPED.1.fq.gz
mv ${sample}_UNMAPPED.2 ${sample}_UNMAPPED.2.fq.gz
gunzip ${sample}_UNMAPPED.2.fq.gz


