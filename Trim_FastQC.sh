#!/bin/bash
#SBATCH --job-name=fastQC_Trim
#SBATCH --mail-user=email
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem-per-cpu=15G
#SBATCH --array=1-17
#SBATCH --output=../Logs/FastQC_Trim_%j.log
#SBATCH --partition=k2-gpu-v100
#SBATCH --chdir=/path/to/file

module load apps/fastqc/0.11.8/noarch
module load apps/trimgalore/0.6.10/noarch

# Specify the path to the config file
config=/path/to/myarray.txt

# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

# Construct the paths to the forward and reverse read files
forward_read="/path/to/files/${sample}_R1_001.fastq.gz"
reverse_read="/path/to/files/${sample}_R2_001.fastq.gz"

echo "This is array task ${SLURM_ARRAY_TASK_ID}, the sample name is ${sample}"

echo ------------------------- FasQC,Trim,FastQC All.fqs --------------------------------

# Provide paths to both forward and reverse read files for paired-end trimming
trim_galore --fastqc -q 30 -length 50 --paired $forward_read $reverse_read