#!/bin/bash
#SBATCH --job-name=Bgla_Genome_Index
#SBATCH --mail-user=email
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem-per-cpu=50G
#SBATCH --output=../Logs/Bgla_Genome_Index
#SBATCH --partition=k2-gpu-v100
#SBATCH --chdir=/path/to/save/files

module load hisat2/2.1.0
module load apps/python/2.7.17/gcc-14.1.0

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/947/242/115/GCF_947242115.1_xgBioGlab47.1/GCF_947242115.1_xgBioGlab47.1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/947/242/115/GCF_947242115.1_xgBioGlab47.1/GCF_947242115.1_xgBioGlab47.1_genomic.gtf.gz

gunzip GCF_947242115.1_xgBioGlab47.1_genomic.fna.gz
gunzip GCF_947242115.1_xgBioGlab47.1_genomic.gtf.gz

mv GCF_947242115.1_xgBioGlab47.1_genomic.fna Bgla_Geno.fa
mv GCF_947242115.1_xgBioGlab47.1_genomic.gtf Bgla_Geno_Annos.gtf

extract_splice_sites.py Bgla_Geno_Annos.gtf >Bgla_Geno.ss
extract_exons.py Bgla_Geno_Annos.gtf >Bgla_Geno.exon
hisat2-build -f -p 4  --ss Bgla_Geno.ss  --exon Bgla_Geno.exon Bgla_Geno.fa Bgla_Geno
