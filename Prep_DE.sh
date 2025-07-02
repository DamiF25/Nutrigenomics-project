#!/bin/bash
#SBATCH --job-name=PrepDE
#SBATCH --mail-user=email
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem-per-cpu=10G
#SBATCH --output=./Logs/PrepDE.log
#SBATCH --partition=bio-compute
#SBATCH --chdir=/path/to/file

module load stringtie/1.3.6
module load apps/python/2.7.8/gcc-4.8.5

./prepDE.py -i ./sample_list.txt -g ./Gene_count_matrix.csv -t ./Transcript_count_matrix.csv -l 130
