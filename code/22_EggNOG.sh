#!/bin/bash -l
#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 00:15:00
#SBATCH -J job_eggNOG
#SBATCH --mail-type=ALL
#SBATCH --mail-user suboticsilvija@yahoo.com
#SBATCH --output=%x.%j.out


# Load modules
module load bioinfo-tools
module load eggNOG-mapper/2.1.9


# Stop executing, if there is error
set -euo pipefail


# Variable
ODir=~/GenomeAnalysis/analysis/22_EggNOG
ANN=~/GenomeAnalysis/analysis/09_Prokka/spades_prokka_noSeq.gff
F=~/GenomeAnalysis/code/sign_genes.fasta


# Commands
emapper.py -i $F --output_dir $ODir --output eggnog -m diamond --cpu 8 --decorate_gff $ANN --override

