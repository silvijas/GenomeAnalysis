#!/bin/bash -l
#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 01:00:00
#SBATCH -J job_synteny
#SBATCH --mail-type=ALL
#SBATCH --mail-user suboticsilvija@yahoo.com
#SBATCH --output=%x.%j.out


# Load modules
module load bioinfo-tools
module load blast/2.15.0+


# Stop executing, if there is error
set -euo pipefail


# Variables
ODir=~/GenomeAnalysis/analysis/08_Synteny
REF=~/GenomeAnalysis/data/REF/GCF_009734005.1_ASM973400v2_genomic.fna.gz
ASM=~/GenomeAnalysis/analysis/04_Spades/contigs.fasta


# Commands
blastn -query $ASM -subject <(gunzip -c $REF) -outfmt 6 -out "${ODir}/spades_blast.out"

