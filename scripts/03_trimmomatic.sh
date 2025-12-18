#!/bin/bash

#SBATCH --job-name=03_trimmomatic
#SBATCH -o 03_trimmomatic_o%j
#SBATCH -e 03_trimmomatic_e%j
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB
#SBATCH --time=72:00:00

source ~/.bash_profile
conda activate rnaseq

helpFunction()
{
   echo ""
   echo "Usage: $0 -F parameterF -R parameterR -K parameterK -S parameterS -L parameterL -T parameterT -C parameterC -H parameterH -M parameterM"
   echo -e "\t-F the file extension for your R1 reads"
   echo -e "\t-R the file extension for your R2 reads"
   echo -e "\t-K parameters for ILLUMINACLIP"
   echo -e "\t-S parameters for SLIDINGWINDOW"
   echo -e "\t-L parameter for LEADING"
   echo -e "\t-T parameter for TRAILING"
   echo -e "\t-C parameter for CROP"
   echo -e "\t-H parameter for HEADCROP"
   echo -e "\t-M parameter for MINLEN"
   exit 1 # Exit script after printing help
}

# ILLUMINACLIP: Cut adapter and other illumina-specific sequences from the read.
#SLIDINGWINDOW: Perform a sliding window trimming, cutting once the average quality within the window falls below a threshold.
#LEADING: Cut bases off the start of a read, if below a threshold quality
#TRAILING: Cut bases off the end of a read, if below a threshold quality
#CROP: Cut the read to a specified length
#HEADCROP: Cut the specified number of bases from the start of the read
#MINLEN: Drop the read if it is below a specified length


while getopts "F:R:K:S:L:T:C:H:M:" opt
do
   case "$opt" in
      F ) parameterF="$OPTARG" ;;
      R ) parameterR="$OPTARG" ;;
      K ) parameterK="$OPTARG" ;;
      S ) parameterS="$OPTARG" ;;
      L ) parameterL="$OPTARG" ;;
      T ) parameterT="$OPTARG" ;;
      C ) parameterC="$OPTARG" ;;
      H ) parameterH="$OPTARG" ;;
      M ) parameterM="$OPTARG" ;;
      #? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$parameterF" ] || [ -z "$parameterR" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

# run command

src=$PWD

mkdir $src/trimmed


for f in $src/raw_data/*$parameterF;
do FBASE=$(basename $f)
	BASE=${FBASE%$parameterF}
	/usr/local/community/Genomics/apps/mambaforge/envs/trimmomatic/bin/trimmomatic PE -threads 4 -phred33 \
	$src/raw_data/${BASE}$parameterF \
	$src/raw_data/${BASE}$parameterR \
	$src/trimmed/${BASE}_trimmed_paired_R1.fastq.gz \
	$src/trimmed/${BASE}_trimmed_unpaired_R1.fastq.gz \
	$src/trimmed/${BASE}_trimmed_paired_R2.fastq.gz \
	$src/trimmed/${BASE}_trimmed_unpaired_R2.fastq.gz \
	$parameterK $parameterS $parameterL $arameterT $parameterC $parameterH $parameterM
done

