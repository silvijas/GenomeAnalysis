#!/bin/bash -l
#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 01:00:00
#SBATCH -J job_spades
#SBATCH --mail-type=ALL
#SBATCH --mail-user suboticsilvija@yahoo.com
#SBATCH --output=%x.%j.out


# Load modules
module load bioinfo-tools
module load spades/3.15.5


# Stop executing, if there is error
set -euo pipefail


# Variables 
ODir=~/GenomeAnalysis/analysis/04_Spades
ILL1=~/GenomeAnalysis/data/DNA/Illumina/E745-1.L500_SZAXPI015146-56_1_clean.fq.gz
ILL2=~/GenomeAnalysis/data/DNA/Illumina/E745-1.L500_SZAXPI015146-56_2_clean.fq.gz


# Commands
zcat ~/GenomeAnalysis/data/DNA/PacBio/* > ~/GenomeAnalysis/total_PacBio.fastq
gzip ~/GenomeAnalysis/total_PacBio.fastq

spades.py -1 $ILL1 -2 $ILL2 --pacbio ~/GenomeAnalysis/total_PacBio.fastq.gz -o $ODir

rm ~/GenomeAnalysis/total_PacBio.fastq.gz
