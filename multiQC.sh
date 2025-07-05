#!/bin/bash
#SBATCH --job-name=MultiQC
#SBATCH --mail-user=email
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem-per-cpu=16G
#SBATCH --output=../Logs/MultiQC.out
#SBATCH --error=../Logs/MultiQC.err
#SBATCH --partition=k2-gpu-v100
#SBATCH --chdir=/path/to/files

module load apps/multiqc/1.15
module load apps/python3/3.12.4/gcc-14.1.0

# Define directories
FASTQC_DIR="/path/to/fastq_files"
MULTIQC_OUT="/path/to/result_output/folder"

# Create output directory if it doesn't exist
mkdir -p "$MULTIQC_OUT"

# Run MultiQC to summarize GC content and other stats
multiqc "$FASTQC_DIR" -o "$MULTIQC_OUT"
