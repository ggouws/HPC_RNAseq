#!/bin/bash

#SBATCH --job-name=kallisto
#SBATCH -o 08_kallisto_o%j
#SBATCH -e 08_kallisto_e%j
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8GB
#SBATCH --time=72:00:00

source ~/.bash_profile
conda activate kallisto


src=$PWD
mkdir $src/kallisto


# Index the transcriptome for pseudoalignment

kallisto index -i $src/reference/transcripts.idx $src/reference/$parameterA*$parameterX

# Quantification

# Loop through each sample and map trimmed paired reads to indexed genome

for f in $src/trimmed/*_trimmed_paired_R1.fastq.gz
do
FBASE=$(basename $f)
BASE=${FBASE%_trimmed_paired_R1.fastq.gz};
kallisto quant \
  -i $src/reference/transcripts.idx \
  -o $src/kallisto/${BASE}_output \
  -b 100 \
  $src/trimmed/${BASE}_trimmed_paired_R1.fastq.gz \
  $src/trimmed/${BASE}_trimmed_paired_R2.fastq.gz \
  -t 4;
mv $src/kallisto/${BASE}_output/abundance.tsv $src/kallisto/${BASE}_abundance.tsv;
mv $src/kallisto/${BASE}_output/run_info.json $src/kallisto/${BASE}_info.json; 
done

rm -r $src/kallisto/*_output
