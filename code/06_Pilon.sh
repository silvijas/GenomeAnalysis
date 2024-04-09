#!/bin/bash -l
#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 01:00:00
#SBATCH -J job_pilon
#SBATCH --mail-type=ALL
#SBATCH --mail-user suboticsilvija@yahoo.com
#SBATCH --output=%x.%j.out


# Load modules
module load bioinfo-tools
module load Pilon/1.24


# Stop executing, if there is error
set -euo pipefail


# Variables
ODir=~/GenomeAnalysis/analysis/06_Pilon
Genome=~/GenomeAnalysis/analysis/03_Flye/assembly.fasta
Frags=~/GenomeAnalysis/analysis/05_Mapping_Illumina2PacBio/aln_Ill2PacBio_sorted.bam

# Commands

## Documentation https://github.com/broadinstitute/pilon/wiki/Requirements-&-Usage
## -- genome : path to the genome sequence that we want to correct.
## -- frags  : path to the bam file containing fragment paired-end alignments, aligned to the --genome argument using bwa.

java -jar $PILON_HOME/pilon.jar --genome $Genome --frags $Frags --output pilon --outdir $ODir 
