#!/bin/bash
#SBATCH --job-name=Parse_Trim_reports
#SBATCH --mail-user=email
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem=512M
#SBATCH --output=../Logs/Parse_Trim_reports.out
#SBATCH --error=../Logs/Parse_Trim_reports.err
#SBATCH --chdir=/path/to/Scripts


# Correct REPORT_DIR
REPORT_DIR="/path/to/Read_Files"
OUTPUT_FILE="parse_trim_sample_summary.tsv"
AVERAGE_FILE="parse_trim_average_summary.tsv"

# Header
echo -e "Sample\tTotal_Reads\tReads_Kept\tBasepairs_Processed\tBasepairs_Kept\tUnpaired_Reads\tPct_Reads_With_Adapters\tPct_Basepairs_Trimmed\tAvg_Read_Length" > $OUTPUT_FILE

# Helper functions
extract_value() {
    grep "$1" "$2" | sed -E 's/[^0-9]*([0-9,]+).*/\1/' | tr -d ','
}
extract_trimmed() {
    grep "$1" "$2" | sed -E 's/[^0-9]*([0-9,]+) bp.*/\1/' | tr -d ','
}
extract_adapter_pct() {
    grep "Reads with adapters" "$1" | sed -E 's/.*\(([0-9.]+)%\).*/\1/'
}

# Process all R1 reports
for R1 in "${REPORT_DIR}"/*_R1_001.fastq.gz_trimming_report.txt; do
    SAMPLE=$(basename "$R1" | sed 's/_R1_001.fastq.gz_trimming_report.txt//')
    R2="${REPORT_DIR}/${SAMPLE}_R2_001.fastq.gz_trimming_report.txt"

    if [[ ! -f "$R2" ]]; then
        echo "R2 report not found for $SAMPLE" >&2
        continue
    fi

    echo "Processing: $SAMPLE"

    # Extract R1 values
    R1_TOTAL=$(extract_value "Total reads processed" "$R1")
    R1_KEPT=$(extract_value "Reads written" "$R1")
    R1_BP_TOTAL=$(extract_value "Total basepairs processed" "$R1")
    R1_BP_KEPT=$(extract_value "Total written" "$R1")
    R1_ADAPT_PCT=$(extract_adapter_pct "$R1")
    R1_BP_TRIMMED=$(extract_trimmed "Quality-trimmed" "$R1")

    # Extract R2 values
    R2_TOTAL=$(extract_value "Total reads processed" "$R2")
    R2_KEPT=$(extract_value "Reads written" "$R2")
    R2_BP_TOTAL=$(extract_value "Total basepairs processed" "$R2")
    R2_BP_KEPT=$(extract_value "Total written" "$R2")
    R2_ADAPT_PCT=$(extract_adapter_pct "$R2")
    R2_BP_TRIMMED=$(extract_trimmed "Quality-trimmed" "$R2")

    # Confirm values
    echo "  R1_TOTAL=$R1_TOTAL, R2_TOTAL=$R2_TOTAL, R1_KEPT=$R1_KEPT, R2_KEPT=$R2_KEPT"

    # Compute metrics with awk
    awk -v sample="$SAMPLE" \
        -v r1_total="$R1_TOTAL" -v r2_total="$R2_TOTAL" \
        -v r1_kept="$R1_KEPT" -v r2_kept="$R2_KEPT" \
        -v r1_bp_total="$R1_BP_TOTAL" -v r2_bp_total="$R2_BP_TOTAL" \
        -v r1_bp_kept="$R1_BP_KEPT" -v r2_bp_kept="$R2_BP_KEPT" \
        -v r1_trim="$R1_BP_TRIMMED" -v r2_trim="$R2_BP_TRIMMED" \
        -v a1="$R1_ADAPT_PCT" -v a2="$R2_ADAPT_PCT" \
        '
   BEGIN {
    total_reads = r1_total + r2_total;
    reads_kept = r1_kept + r2_kept;
    bp_total = r1_bp_total + r2_bp_total;
    bp_kept = r1_bp_kept + r2_bp_kept;
    unpaired = total_reads - 2 * r1_kept;
    pct_adapt = (a1 + a2) / 2;
    pct_trim = ((r1_trim + r2_trim) / bp_total) * 100;
    avg_read_length = bp_total / total_reads;
    printf "%s\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.1f\n", sample, total_reads, reads_kept, bp_total, bp_kept, unpaired, pct_adapt, pct_trim, avg_read_length;
}
done

# Compute overall averages
tail -n +2 "$OUTPUT_FILE" | awk -F'\t' '
{
    total_reads += $2
    reads_kept += $3
    bp_total += $4
    bp_kept += $5
    unpaired += $6
    pct_adapt += $7
    pct_trimmed += $8
    avg_read_length += $9
    count += 1
}
END {
    if (count > 0)
        printf "AVERAGE\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.2f\t%.2f\t%.1f\n", total_reads/count, reads_kept/count, bp_total/count, bp_kept/count, unpaired/count, pct_adapt/count, pct_trimmed/count, avg_read_length/count
    else
        print "No data rows to summarize"
}' >> "$AVERAGE_FILE"

echo "Summary written to $OUTPUT_FILE"
echo "Averages written to $AVERAGE_FILE"