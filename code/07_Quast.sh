#!/bin/bash -l
#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 01:00:00
#SBATCH -J job_quast
#SBATCH --mail-type=ALL
#SBATCH --mail-user suboticsilvija@yahoo.com
#SBATCH --output=%x.%j.out


# Load modules
module load bioinfo-tools
module load quast/5.0.2


# Stop executing, if there is error
set -euo pipefail


# Variables
ODir=~/GenomeAnalysis/analysis/07_Quast
REF=~/GenomeAnalysis/data/REF/GCF_009734005.1_ASM973400v2_genomic.fna.gz

## Spades assembly
S_ASS=~/GenomeAnalysis/analysis/04_Spades/contigs.fasta
P_ASS=~/GenomeAnalysis/analysis/06_Pilon/pilon.fasta

# Commands 
quast.py $S_ASS -r $REF -o "${ODir}/Spades" --threads 8 --gene-finding
quast.py $P_ASS -r $REF -o "${ODir}/Pilon"  --threads 8 --gene-finding
