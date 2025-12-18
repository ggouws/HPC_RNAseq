#!/bin/bash

#SBATCH --job-name=07_htseq-count
#SBATCH -o 07_htseq-count_o%j
#SBATCH -e 07_htseq-count_e%j
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8GB
#SBATCH --time=72:00:00

source ~/.bash_profile
conda activate rnaseq

helpFunction()
{
  echo ""
  echo "Usage: $0 -m parameterM -s parameterS -t parameterT -i parameterI -r parameterR -A parameterA -X parameterX"
  echo -e "\t-m mode for handling reads overlapping more than one genomic feature"
  echo -e "\t-s indication of whether the data are stranded ('yes') or not ('no')
  echo -e "\t-t the feature to be used (e.g., gene or exon)"
  echo -e "\t-i GFF attribute used to identify the feature"
  echo -e "\t-r alignment sorted by read name ('name') or alignment position ('pos')"
  echo -e "\t-A start of the name or accession number of genomic/transcriptomic annotation file"
  echo -e "\t-X extension of the the genomic/transcriptomic annotation file"
  exit 1 # Exit script after printing help
}

while getopts "m:s:t:i:r:A:X:" opt
do
   case "$opt" in
      m ) parameterM="$OPTARG" ;;
      s ) parameterS="$OPTARG" ;;
      t ) parameterT="$OPTARG" ;;
      i ) parameterI="$OPTARG" ;;
      r ) parameterR="$OPTARG" ;;
      A ) parameterA="$OPTARG" ;;
      X ) parameterX="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$parameterM" ] || [ -z "$parameterS" ] || [ -z "$parameterT" ] || [ -z "$parameterI" ] || [ -z "$parameterR" ] || [ -z "$parameterA" ] || [ -z "$parameterX" ]
then
   echo "One or more of the mandatory arguments is/are missing. Please check these your command line submission.";
   helpFunction
fi

src=$PWD
mkdir $src/htseq


# Loop through each sample and map trimmed paired reads to indexed genome

for f in $src/aligned_clean/*.clean.bam
do
FBASE=$(basename $f)
BASE=${FBASE%.clean.bam};
htseq-count \
-m "$parameterM" \
-s "$parameterS" \
-t "$parameterT" \
-i "$parameterI" \
-r "$parameterR" \
-f bam \
"$src/aligned_clean/${BASE}.clean.bam" "$src/reference/"$parameterA"*"$parameterX" > \
"$src/htseq/${BASE}_"$parameterT"count.htseq";
done
