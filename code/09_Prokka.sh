#!/bin/bash -l
#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 01:00:00
#SBATCH -J job_prokka
#SBATCH --mail-type=ALL
#SBATCH --mail-user suboticsilvija@yahoo.com
#SBATCH --output=%x.%j.out


# Load modules
module load bioinfo-tools
module load prokka/1.45-5b58020


# Stop executing, if there is error
set -euo pipefail


# Varibles
ODir=~/GenomeAnalysis/analysis/09_Prokka
ASM=~/GenomeAnalysis/analysis/04_Spades/contigs.fasta


# Commands
prokka --outdir $ODir --force --prefix spades_prokka \
       --genus Enterococcus --species faecium --strain E745 \
       --kingdom Bacteria --gram pos \
       $ASM
            
