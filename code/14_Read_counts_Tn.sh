#!/bin/bash -l
#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 09:00:00
#SBATCH -J job_read_counts_Tn
#SBATCH --mail-type=ALL
#SBATCH --mail-user suboticsilvija@yahoo.com
#SBATCH --output=%x.%j.out


# Load modules
module load bioinfo-tools
module load htseq/2.0.2


# Stop executing, if there is error
set -euo pipefail


# Variables
ODir=~/GenomeAnalysis/analysis/14_Read_counts_Tn_3_eg
BAMs=~/GenomeAnalysis/analysis/13_Tn_mapping
#ANN=~/GenomeAnalysis/analysis/09_Prokka/spades_prokka_noSeq_90perc.gff
ANN=~/GenomeAnalysis/analysis/22_EggNOG/eggnog.emapper.decorated.gff

# Commands

## Function that assigns tag to file depending on a folder where it originated.
get_folder_tag(){
	if [[ $1 == *"BHI"* ]]; then
		echo "BHI"
	elif [[ $1 == *"HSerum"* ]]; then 
		echo "HSerum"           
        else
		echo "Serum"
	fi
}

for file in ${BAMs}/*_sorted.bam; do

	echo $file

	id=$(echo $file | grep -oP 'ERR[0-9]+')
	echo $id

	folder_tag=$(get_folder_tag "$file")
	echo $folder_tag

	htseq-count -t CDS -r pos -i ID -f bam -s no $file $ANN > "${ODir}/${folder_tag}_${id}_readcount.txt"
done




