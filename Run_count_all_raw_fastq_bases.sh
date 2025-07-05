#!/bin/bash
#SBATCH --job-name=Count_Rwa_FASTQ_Bases
#SBATCH --mail-user=email
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem-per-cpu=4G
#SBATCH --output=../Logs/Count_Raw_FASTQ_Bases.out
#SBATCH --error=../Logs/Count_Raw_FASTQ_Bases.err
#SBATCH --partition=k2-gpu-v100
#SBATCH --chdir=/path/to/files

# Go to the directory containing the FASTQ files
cd /path/to/files || exit 1

# Create output CSV
output_file="merged_raw_fastq_base_counts.csv"
echo "Sample,Total_Bases,Raw_Data_Gb" > "$output_file"

# Loop through all R1 files
for r1 in *_R1_001.fastq.gz; do
    sample=$(basename "$r1" _R1_001.fastq.gz)
    r2="${sample}_R2_001.fastq.gz"
  
 if [[ -f "$r2" ]]; then
        echo "Processing $sample..."

        bases_r1=$(zcat "$r1" | awk 'NR % 4 == 2 {bases += length($0)} END {print bases}')
        bases_r2=$(zcat "$r2" | awk 'NR % 4 == 2 {bases += length($0)} END {print bases}')
        total_bases=$((bases_r1 + bases_r2))
        gb=$(awk -v b="$total_bases" 'BEGIN {printf "%.3f", b / 1e9}')

        echo "$sample,$total_bases,$gb" >> "$output_file"
    else
        echo ⚠️ Warning: R2 file missing for $sample"
    fi
done

echo "✅ Done. Output written to $output_file"
