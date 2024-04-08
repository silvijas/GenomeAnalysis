#!/bin/bash -l
#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 03:00:00
#SBATCH -J job_trimmomatic
#SBATCH --mail-type=ALL
#SBATCH --mail-user suboticsilvija@yahoo.com
#SBATCH --output=%x.%j.out


# Load modules
#module load bioinfo-tools
#module load trimmomatic/0.39
#module load java/sun_jdk1.8.0_151


# Commands
java -jar /sw/bioinfo/trimmomatic/0.39/snowy/trimmomatic-0.39.jar PE -threads 2 -phred33 ../data/DNA/Illumina/E745-1.L500_SZAXPI015146-56_1_clean.fq.gz ../data/DNA/Illumina/E745-1.L500_SZAXPI015146-56_2_clean.fq.gz ../results/01_Trimmomatic/E745-1_out_forward_paired.fastq ../results/01_Trimmomatic/E745-1_out_forward_unpaired.fastq ../results/01_Trimmomatic/E745-1_out_reverse_paired.fastq ../results/01_Trimmomatic/E745-1_out_reverse_unpaired.fastq ILLUMINACLIP:/sw/bioinfo/trimmomatic/0.39/snowy/adapters/TruSeq3-PE.fa:2:30:7 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:36

