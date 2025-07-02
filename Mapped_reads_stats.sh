#!/bin/bash
#SBATCH --job-name=mapped_reads_stats
#SBATCH --mail-user=email
#SBATCH --output=logs/mapStats_%A_%a.out
#SBATCH --error=logs/mapStats_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-17 
#SBATCH --partition=k2-gpu-v100
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --chdir=/path/to/result_files

module load apps/samtools/1.9/gcc-14.1.0

# Define input directory and output
INPUT_DIR=/file/to/files
OUTPUT_FILE="mapped_unmapped_summary.tsv"

# Make sure logs directory exists
mkdir -p logs

# Get list of SAM files
SAM_FILES=(${INPUT_DIR}/*.sam)
SAM_FILE=${SAM_FILES[$SLURM_ARRAY_TASK_ID]}
BASENAME=$(basename "$SAM_FILE" .sam)

# Convert SAM to BAM and keep only primary alignments (removes secondary/supplementary alignments)
samtools view -b -F 0x100 -F 0x800 "$SAM_FILE" > "${BASENAME}.primary.bam"

# Run stats only on primary alignments
STATS=$(samtools flagstat "${BASENAME}.primary.bam")

# Extract relevant numbers
TOTAL=$(echo "$STATS" | grep "in total" | awk '{print $1}')
MAPPED=$(echo "$STATS" | grep "mapped (" | awk '{print $1}')
UNMAPPED=$((TOTAL - MAPPED))
PCT_MAPPED=$(awk -v a=$MAPPED -v b=$TOTAL 'BEGIN {printf "%.2f", (a/b)*100}')
PCT_UNMAPPED=$(awk -v a=$UNMAPPED -v b=$TOTAL 'BEGIN {printf "%.2f", (a/b)*100}')

# Print header only once
if [ "$SLURM_ARRAY_TASK_ID" -eq 1 ]; then
    echo -e "Sample\tTotal_Reads\tMapped\tUnmapped\tPct_Mapped\tPct_Unmapped" > "$OUTPUT_FILE"
fi

# Output results
echo -e "${BASENAME}\t${TOTAL}\t${MAPPED}\t${UNMAPPED}\t${PCT_MAPPED}\t${PCT_UNMAPPED}" >> "$OUTPUT_FILE"

# Clean up intermediate BAM
rm "${BASENAME}.bam"