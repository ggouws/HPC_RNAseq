#!/bin/bash

#SBATCH --job-name=01_reference_download
#SBATCH -o 01_reference_download.sh_o%j
#SBATCH -e 01_reference_download.sh_e%j
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB
#SBATCH --time=00:20:00

source ~/.bash_profile
conda activate rnaseq

src=$PWD

helpFunction()
{
   echo ""
   echo "Usage: $0 -L parameterL"
   echo -e "\t-L the NCBI-datasets download link for your genome/transcriptome"
   exit 1 # Exit script after printing help
}

while getopts "L:" opt
do
   case "$opt" in
      L ) parameterL="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$parameterL" ]
then
   echo "A NCBI-datasets link was not provided";
   helpFunction
fi

# Begin script in case all parameters are correct

mkdir $src/reference
cd $src/reference

### Download genome/transcriptome
eval $parameterL

# Unzip the downloaded data, move genome/transcriptome and annotation files and tidy (remove downloaded zip)
unzip ncbi_dataset
cp ncbi_dataset/data/*/*.fna ncbi_dataset/data/*/*.fasta ncbi_dataset/data/*/*.gff ncbi_dataset/data/*/*.gff3 ncbi_dataset/*/*.gtf .
rm -r ncbi_*
rm md5sum.txt
rm README.md
