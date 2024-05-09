#!/bin/bash -l
#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 01:00:00
#SBATCH -J job_flye
#SBATCH --mail-type=ALL
#SBATCH --mail-user suboticsilvija@yahoo.com
#SBATCH --output=%x.%j.out


# Load modules
module load bioinfo-tools
module load Flye/2.9.1


# Variables
DATA=~/GenomeAnalysis/data/DNA/PacBio/*
ODir=~/GenomeAnalysis/analysis/03_Flye/


# Commands
flye --pacbio-raw $DATA --out-dir $ODir --threads 4
