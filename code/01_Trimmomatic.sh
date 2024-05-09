#!/bin/bash -l
#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 03:00:00
#SBATCH -J job_trimmomatic
#SBATCH --mail-type=ALL
#SBATCH --mail-user suboticsilvija@yahoo.com
#SBATCH --output=%x.%j.out


# Load modules
module load bioinfo-tools
module load trimmomatic/0.39
module load java/sun_jdk1.8.0_151


# Variables 
ILL1=~/GenomeAnalysis/data/DNA/Illumina/E745-1.L500_SZAXPI015146-56_1_clean.fq.gz
ILL2=~/GenomeAnalysis/data/DNA/Illumina/E745-1.L500_SZAXPI015146-56_2_clean.fq.gz

ODir=~/GenomeAnalysis/analysis/01_Trimmomatic
FP=${ODir}/E745-1_out_forward_paired.fastq
RP=${ODir}/E745-1_out_reverse_paired.fastq

FUP=${ODir}/E745-1_out_forward_unpaired.fastq
RUP=${ODir}/E745-1_out_reverse_unpaired.fastq

LOG=${ODir}/log.txt


# Commands
java -jar /sw/bioinfo/trimmomatic/0.39/snowy/trimmomatic-0.39.jar PE -threads 2 -phred33 $ILL1 $ILL2 $FP $FUP $RP $RUP ILLUMINACLIP:/sw/bioinfo/trimmomatic/0.39/snowy/adapters/TruSeq3-PE.fa:2:30:7 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:36 > $LOG

