#!/bin/bash

#SBATCH --job-name=07_htseq-count
#SBATCH -o 07_htseq-count_o%j
#SBATCH -e 07_htseq-count_e%j
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH -A molecolb
#SBATCH -p molecolb
#SBATCH --mem-per-cpu=8GB
#SBATCH --time=72:00:00

source ~/.bash_profile
#####GET THE HTSEQ ENVIRONMENT
conda activate hisat2

helpFunction()
{
  echo ""
  echo "Usage: $0 -m parameterM -t parameterT -i parameterI -r parameterR"
  echo -e "\t-m mode for handling reads overlapping more than one genomic feature"
  echo -e "\t-t the feature to be used (e.g., gene or exon)"
  echo -e "\t-i GFF attribute used to identify the feature"
  echo -e "\t-r alignment sorted by read name ('name') or alignment position ('pos')
  exit 1 # Exit script after printing help
}

while getopts "m:t:i:r:" opt
do
   case "$opt" in
      m ) parameterA="$OPTARG" ;;
      t ) parameterT="$OPTARG" ;;
      i ) parameterI="$OPTARG" ;;
      r ) parameterR="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$parameterA" ] || [ -z "$parameterT" ] || [ -z "$parameterI" ] || [ -z "$parameterR" ]
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
-m parameterM \
-s yes \
-t parameterT \
-i parameterI \
-r parameterR \
-f bam \
$src/aligned_clean/${BASE}.clean.bam \
$src/reference/*GFF_FILE_HERE > \
$src/htseq/${BASE}_genecount.htseq;
done
