#!/bin/bash -l
#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 01:00:00
#SBATCH -J job_FastQC_Illumina
#SBATCH --mail-type=ALL
#SBATCH --mail-user suboticsilvija@yahoo.com
#SBATCH --output=%x.%j.out


# Load modules
module load bioinfo-tools
module load FastQC/0.11.9 


# Variables
ILL1=~/GenomeAnalysis/data/DNA/Illumina/E745-1.L500_SZAXPI015146-56_1_clean.fq.gz
ILL2=~/GenomeAnalysis/data/DNA/Illumina/E745-1.L500_SZAXPI015146-56_2_clean.fq.gz


# Commands
fastqc -o $ILL1 $ILL2 

