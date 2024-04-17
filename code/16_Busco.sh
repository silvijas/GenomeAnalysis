#!/bin/bash -l
#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 02:00:00
#SBATCH -J job_busco_pilon
#SBATCH --mail-type=ALL
#SBATCH --mail-user suboticsilvija@yahoo.com
#SBATCH --output=%x.%j.out


# Load modules
module load bioinfo-tools
module load BUSCO/5.5.0


# Stop executing, if there is error
set -euo pipefail


# Variables
ODir=~/GenomeAnalysis/analysis/16_Busco
REF=~/GenomeAnalysis/data/REF/GCF_009734005.1_ASM973400v2_genomic.fna.gz
REF_f=~/GenomeAnalysis/data/REF/GCF_009734005.1_ASM973400v2_genomic.fna

## Spades assembly
S_ASS=~/GenomeAnalysis/analysis/04_Spades/contigs.fasta
P_ASS=~/GenomeAnalysis/analysis/06_Pilon/pilon.fasta


# Commands 

source $AUGUSTUS_CONFIG_COPY

busco -i $S_ASS -o spades_bacteria        -l $BUSCO_LINEAGE_SETS/bacteria               -m genome -c 8 -f --out_path $ODir --download_path $ODir
busco -i $P_ASS -o pilon_bacteria         -l $BUSCO_LINEAGE_SETS/bacteria               -m genome -c 8 -f --out_path $ODir --download_path $ODir

busco -i $S_ASS -o spades_e_bacterales    -l $BUSCO_LINEAGE_SETS/enterobacterales_odb10 -m genome -c 8 -f --out_path $ODir --download_path $ODir
busco -i $P_ASS -o pilon_e_bacterales     -l $BUSCO_LINEAGE_SETS/enterobacterales_odb10 -m genome -c 8 -f --out_path $ODir --download_path $ODir

busco -i $REF_f -o reference_e_bacterales -l $BUSCO_LINEAGE_SETS/enterobacterales_odb10 -m genome -c 8 -f --out_path $ODir --download_path $ODir

