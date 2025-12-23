#!/bin/bash

#SBATCH --job-name=05_reference_mapping
#SBATCH -o 05_reference_mapping_o%j
#SBATCH -e 05_reference_mapping_e%j
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4GB
#SBATCH --time=72:00:00

source ~/.bash_profile
conda activate rnaseq

helpFunction()
{
  echo ""
  echo "Usage: $0 -A parameterA -X parameterX"
  echo -e "\t-A the accession number of your reference genome or transcriptome"
  echo -e "\t-X the file extension of your reference genome or transcriptome"
  exit 1 # Exit script after printing help
}

while getopts "A:X:" opt
do
   case "$opt" in
      A ) parameterA="$OPTARG" ;;
      X ) parameterX="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$parameterA" ] || [ -z "$parameterX" ]
then
   echo "You have not provided the accession or/and extension of your reference genome/transcriptome file. Please check these.";
   helpFunction
fi

src=$PWD
mkdir $src/aligned

# Index the reference genome

hisat2-build -f $src/reference/$parameterA*$parameterX $src/reference/idxref

# Loop through each sample and map trimmed paired reads to indexed genome

for f in $src/trimmed/*_trimmed_paired_R1.fastq.gz
do
FBASE=$(basename $f)
BASE=${FBASE%_trimmed_paired_R1.fastq.gz};
hisat2 \
  -x $src/reference/idxref \
  -1 $src/trimmed/${BASE}_trimmed_paired_R1.fastq.gz \
  -2 $src/trimmed/${BASE}_trimmed_paired_R2.fastq.gz \
  -p 8 \
  -q \
  --met-file $src/aligned/${BASE}.stats | \
  samtools sort -O BAM > $src/aligned/${BASE}.bam;
echo ${BASE} >> $src/quality_reports/initial_mapping_quality;
echo "" >> $src/quality_reports/initial_mapping_quality;
samtools flagstat $src/aligned/${BASE}.bam >> $src/quality_reports/initial_mapping_quality;
echo "" >> $src/quality_reports/initial_mapping_quality;
echo "" >> $src/quality_reports/initial_mapping_quality;
done
