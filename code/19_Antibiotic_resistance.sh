#!/bin/bash -l
#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 02:00:00
#SBATCH -J job_contig_blast
#SBATCH --mail-type=ALL
#SBATCH --mail-user suboticsilvija@yahoo.com
#SBATCH --output=%x.%j.out


# Load modules
module load bioinfo-tools
module load abricate/1.0.1-20200512-955d402 blast/2.15.0+ blast_databases


# Stop executing, if there is error
set -euo pipefail


# Variables 
ODir=~/GenomeAnalysis/analysis/19_Antibiotic_resistance
ASMB=~/GenomeAnalysis/analysis/04_Spades/contigs.fasta


# Commands

abricate --db ncbi $ASMB > "${ODir}/abricate_output.tsv"

## From the identified sequences contig 4 was not identified as plasmid in the plasmid identification.
##   Blast was performed to try to identify it is chromosome or plasmid.
blastn -db nt -query "${ODir}/contig_04.fasta" > "${ODir}/contig_04_blast.out"

