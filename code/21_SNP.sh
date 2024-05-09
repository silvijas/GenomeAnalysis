#!/bin/bash -l
#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 08:00:00
#SBATCH -J job_SNP
#SBATCH --mail-type=ALL
#SBATCH --mail-user suboticsilvija@yahoo.com
#SBATCH --output=%x.%j.out



# Stop executing, if there is error
set -euo pipefail


# Load modules
module load bioinfo-tools
module load python/3.9.5 texlive/2022-09-18
module load bwa/0.7.17 samtools picard/3.1.1 bcftools 


# Variables
ODir=~/GenomeAnalysis/analysis/21_SNP
FQs=~/GenomeAnalysis/data/SNP
REF=~/GenomeAnalysis/analysis/04_Spades/contigs.fasta
E=.fastq.gz

# Commands
declare -A ids

for file in ${FQs}/*.fastq.gz; do
    	id=$(echo $file | grep -oP 'SRR[0-9]+')
#    	echo "Extracted ID: $id from $file"

	ids["$id"]=1
done

# Output file
OFile="${ODir}/SNPs_per_sample.tsv" > "$OFile"

for id in "${!ids[@]}"; do
	echo $id    

        ILL1="${FQs}/${id}_1${E}"
        ILL2="${FQs}/${id}_2${E}"

	echo $ILL1
	echo $ILL2

        prefix=${ODir}/${id}
        echo $prefix

        bwa mem -t 8 $REF $ILL1 $ILL2 > ${prefix}.sam

        samtools view  -b ${prefix}.sam | samtools sort -o  ${prefix}_sort.bam -
	rm ${prefix}.sam

	java -jar $PICARD MarkDuplicates I=${prefix}_sort.bam O=${prefix}_sort_dupl.bam M=${prefix}_dupl_metrics.txt REMOVE_DUPLICATES=false	

 	samtools index ${prefix}_sort_dupl.bam -o  ${prefix}_sort_dupl.bai

        
	bcftools mpileup -Ou -f $REF ${prefix}_sort_dupl.bam | bcftools call -mv -Ov -o ${prefix}_variants.vcf


	bcftools filter -s LOWQUAL -e 'QUAL<30 || DP<10 || MQ<30' -m '+' -O v -o ${prefix}_variants_filtered.vcf ${prefix}_variants.vcf


	count=$(awk '$7 == "PASS"' ${prefix}_variants_filtered.vcf | wc -l)
	echo -e "$id\t$count" >> "$OFile"

done
