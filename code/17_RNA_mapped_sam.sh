#!/bin/bash -l
#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 09:00:00
#SBATCH -J job_read_counts
#SBATCH --mail-type=ALL
#SBATCH --mail-user suboticsilvija@yahoo.com
#SBATCH --output=%x.%j.out


# Load modules
module load bioinfo-tools
module load samtools/1.19


# Stop executing, if there is error
set -euo pipefail


# Variables
ODir=~/GenomeAnalysis/analysis/17_RNA_mapped_sam
BAMs=~/GenomeAnalysis/analysis/10_RNA_mapping
ANN=~/GenomeAnalysis/analysis/09_Prokka/spades_prokka_noSeq.gff

# Commands

## Function that assigns tag to file depending on a folder where it originated.
get_folder_tag(){
	if [[ $1 == *"BH"* ]]; then
		echo "BH"
	elif [[ $1 == *"Serum"* ]]; then 
		echo "Serum"
	else 
		echo "Error"
	fi
}

for file in ${BAMs}/*_sorted.bam; do

#	echo $file

	id=$(echo $file | grep -oP 'ERR[0-9]+')
	echo $id

	folder_tag=$(get_folder_tag "$file")

	samtools flagstat $file -@ 8 -O tsv > "${ODir}/${folder_tag}_${id}_stats.tsv"
done




