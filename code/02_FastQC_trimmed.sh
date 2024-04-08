#!/bin/bash -l
#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 01:00:00
#SBATCH -J job_
#SBATCH --mail-type=ALL
#SBATCH --mail-user suboticsilvija@yahoo.com
#SBATCH --output=%x.%j.out


# Load modules
module load bioinfo-tools
module load FastQC/0.11.9


# Your commands
fastqc -o ~/GA/GenomeAnalysis/results/02_FastQC_trimmed/ ~/GA/GenomeAnalysis/results/01_Trimmomatic/E745-1_out_forward_paired.fastq ~/GA/GenomeAnalysis/results/01_Trimmomatic/E745-1_out_forward_unpaired.fastq ~/GA/GenomeAnalysis/results/01_Trimmomatic/E745-1_out_reverse_paired.fastq ~/GA/GenomeAnalysis/results/01_Trimmomatic/E745-1_out_reverse_unpaired.fastq
