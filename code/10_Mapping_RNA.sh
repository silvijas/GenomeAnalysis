#!/bin/bash -l
#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 02:00:00
#SBATCH -J job_RNAmapping
#SBATCH --mail-type=ALL
#SBATCH --mail-user suboticsilvija@yahoo.com
#SBATCH --output=%x.%j.out


# Load modules
module load bioinfo-tools
module load bwa/0.7.17 samtools/1.19


# Stop executing, if there is error
set -euo pipefail


# Variables
ODir=~/GenomeAnalysis/analysis/10_RNA_mapping

ASMB=~/GenomeAnalysis/analysis/04_Spades/contigs.fasta

FQs_B=~/GenomeAnalysis/data/RNA/RNA-Seq_BH
FQs_S=~/GenomeAnalysis/data/RNA/RNA-Seq_Serum


# Commands
bwa index $ASMB
 

## Arrays
declare -A ids		# It will contain IDs of samples.
declare -A dict		# It will contain IDs, Serum/BH, directory of sample.


## Function that assigns tag to file depending on a folder where it originated.
get_folder_tag(){
	if [[ $1 == *"Seq_BH"* ]]; then
		echo "BH"
	elif [[ $1 == *"Seq_Serum"* ]]; then 
		echo "Serum"
	else 
		echo "Error"
	fi
}


for dir in "$FQs_B" "$FQs_S"; do	
	for file in "${dir}"/*paired*.fastq.gz; do
#		echo $file
 
		id=$(echo $file | grep -oP 'ERR[0-9]+')		
		ids[$id]=1

		folder_tag=$(get_folder_tag "$file")
                
		dict["$id,$folder_tag,$dir"]=1
	done
done



ST=/trim_paired_
E1=_pass_1.fastq.gz
E2=_pass_2.fastq.gz

echo "sample_id,count_mapped_reads,count_total_reads,percentage" > "${ODir}/percentages.csv"

for key in "${!dict[@]}"; do
	IFS=',' read -r id folder_tag dir <<< "$key"
	echo "id= $id, folder_tag= $folder_tag, folder= $dir"

	ILL1="${dir}${ST}${id}${E1}"
	ILL2="${dir}${ST}${id}${E2}"

	prefix=${ODir}/${folder_tag}_${id}_paired
	echo $prefix

	bwa mem -t 8 $ASMB $ILL1 $ILL2 > ${prefix}.sam	
	samtools view  ${prefix}.sam 	    -bo ${prefix}.bam
	samtools sort  ${prefix}.bam 	    -o 	${prefix}_sorted.bam
	samtools index ${prefix}_sorted.bam -o  ${prefix}_sorted.bai

	rm ${prefix}.sam
 	rm ${prefix}.bam


	# For assesing quality of mapping
	samtools flagstat ${prefix}_sorted.bam > ${prefix}_flagstat.txt
	samtools idxstats ${prefix}_sorted.bam > ${prefix}_idxstats.txt

	total_reads=$(grep -P 'in total' ${prefix}_flagstat.txt | grep -oP '\d+' | head -1)
    	mapped_reads=$(grep -P 'mapped \(' ${prefix}_flagstat.txt | grep -oP '\d+' | head -1)
 	percentage=$(grep -P 'mapped \(' ${prefix}_flagstat.txt | head -1 | grep -oP '(?<=\()\d+\.\d+(?=%)')

    	echo "$id,$mapped_reads,$total_reads,$percentage" >> "${ODir}/percentages.csv"	
done



 
