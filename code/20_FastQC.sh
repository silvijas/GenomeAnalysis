#!/bin/bash -l
#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 08:00:00
#SBATCH -J job_FastQC
#SBATCH --mail-type=ALL
#SBATCH --mail-user suboticsilvija@yahoo.com
#SBATCH --output=%x.%j.out



# Stop executing, if there is error
set -euo pipefail


# Load modules
module load bioinfo-tools
module load FastQC/0.11.9


# Variables
ODir=~/GenomeAnalysis/analysis/20_FastQC
FQs=~/GenomeAnalysis/data/SNP
REF=~/GenomeAnalysis/analysis/04_Spades/contigs.fasta


# Commands
declare -A ids

E=.fastq.gz

for file in ${FQs}/*.fastq.gz; do
    	id=$(echo $file | grep -oP 'SRR[0-9]+')
#    	echo "Extracted ID: $id from $file"

	ids["$id"]=1
done

for id in "${!ids[@]}"; do
	echo $id    

	ODir_id="${ODir}/${id}"
	echo $ODir_id	

        ILL1="${FQs}/${id}_1${E}"
        ILL2="${FQs}/${id}_2${E}"

        echo $ILL1
        echo $ILL2
	
	fastqc $ILL1 $ILL2 -o $ODir	
done
