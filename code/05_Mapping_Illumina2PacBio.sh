#!/bin/bash -l
#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 02:00:00
#SBATCH -J job_03_convert_sort_index
#SBATCH --mail-type=ALL
#SBATCH --mail-user suboticsilvija@yahoo.com
#SBATCH --output=%x.%j.out


# Load modules
module load bioinfo-tools
module load bwa/0.7.17 bwa-mem2 samtools/1.19


# Stop executing, if there is error
set -euo pipefail


# Variables
ODir=~/GenomeAnalysis/analysis/05_Mapping_Illumina2PacBio

ASMB=~/GenomeAnalysis/analysis/03_Flye/assembly.fasta
ILL1=~/GenomeAnalysis/data/DNA/Illumina/E745-1.L500_SZAXPI015146-56_1_clean.fq.gz
ILL2=~/GenomeAnalysis/data/DNA/Illumina/E745-1.L500_SZAXPI015146-56_2_clean.fq.gz


# Commands

bwa index $ASMB 

bwa mem -t 8 $ASMB $ILL1 $ILL2 > "${ODir}/aln_Ill2PacBio.sam"

## Converting file SAM -> BAM	From :: https://www.htslib.org/doc/samtools-view.html
samtools view -bo "${ODir}/aln_Ill2PacBio.bam" "${ODir}/aln_Ill2PacBio.sam"

## Sorting and indexing sorted BAM file, as it is generally the preferred format for further analysis.
samtools sort "${ODir}/aln_Ill2PacBio.bam" -o "${ODir}/aln_Ill2PacBio_sorted.bam"
samtools index "${ODir}/aln_Ill2PacBio_sorted.bam"
