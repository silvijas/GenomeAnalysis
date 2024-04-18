#!/bin/bash -l
#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 02:00:00
#SBATCH -J job_resistance
#SBATCH --mail-type=ALL
#SBATCH --mail-user suboticsilvija@yahoo.com
#SBATCH --output=%x.%j.out


# Load modules
module load bioinfo-tools
module load abricate/1.0.1-20200512-955d402


# Stop executing, if there is error
set -euo pipefail


# Variables 
ODir=~/GenomeAnalysis/analysis/19_Antibiotic_resistance
ASMB=~/GenomeAnalysis/analysis/04_Spades/contigs.fasta


# Commands
abricate --db ncbi $ASMB > "${ODir}/abricate_output.tsv"
