#!/bin/bash -l
#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 03:00:00
#SBATCH -J job_TNmapping
#SBATCH --mail-type=ALL
#SBATCH --mail-user suboticsilvija@yahoo.com
#SBATCH --output=%x.%j.out


# Load modules
module load bioinfo-tools
module load bwa/0.7.17 samtools/1.19


# Stop executing, if there is error
set -euo pipefail


# Variables
ODir=~/GenomeAnalysis/analysis/13_Tn_mapping

ASMB=~/GenomeAnalysis/analysis/04_Spades/contigs.fasta

FQs_B=~/GenomeAnalysis/data/RNA/Tn-Seq_BHI
FQs_S=~/GenomeAnalysis/data/RNA/Tn-Seq_Serum
FQs_H=~/GenomeAnalysis/data/RNA/Tn-Seq_HSerum

# Commands

## Already indexed for RNA mapping
bwa index $ASMB
 

## Arrays
declare -A ids		# It will contain IDs of samples.
declare -A dict		# It will contain IDs, Serum/BH, directory of sample.


## Function that assigns tag to file depending on a folder where it originated.
get_folder_tag(){
	if [[ $1 == *"Seq_BH"* ]]; then
		echo "BHI"
	elif [[ $1 == *"Seq_Serum"* ]]; then 
		echo "Serum"
	elif [[ $1 == *"Seq_HSerum"* ]]; then
		echo "HSerum"
	else 
		echo "Error"
	fi
}


for dir in "$FQs_B" "$FQs_S" "$FQs_H"; do	
	for file in "${dir}"/*.fastq.gz; do
#		echo $file
 
		id=$(echo $file | grep -oP 'ERR[0-9]+')		
		ids[$id]=1

		folder_tag=$(get_folder_tag "$file")
                
		dict["$id,$folder_tag,$dir"]=1
	done
done



ST=/trim_
E1=_pass.fastq.gz

for key in "${!dict[@]}"; do
	IFS=',' read -r id folder_tag dir <<< "$key"
#	echo "id= $id, folder_tag= $folder_tag, folder= $dir"

	ILL1="${dir}${ST}${id}${E1}"

	prefix=${ODir}/${folder_tag}_${id}
#	echo $prefix

	bwa mem -t 8 $ASMB $ILL1 > ${prefix}.sam	
	samtools view  ${prefix}.sam 	    -bo ${prefix}.bam
	samtools sort  ${prefix}.bam 	    -o 	${prefix}_sorted.bam
	samtools index ${prefix}_sorted.bam -o  ${prefix}_sorted.bai

	rm ${prefix}.sam
	rm ${prefix}.bam
done


 

