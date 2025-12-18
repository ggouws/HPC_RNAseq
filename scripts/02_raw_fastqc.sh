#!/bin/bash

#SBATCH --job-name=02_raw_fastqc
#SBATCH -o 02_raw_fastqc_o%j
#SBATCH -e 02_raw_fastqc_e%j
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2GB
#SBATCH --time=48:00:00

source ~/.bash_profile
conda activate rnaseq

helpFunction()
{
   echo ""
   echo "Usage: $0 -F parameterF -R parameterR"
   echo -e "\t-F the file extension for your R1 reads"
   echo -e "\t-R the file extension for your R2 reads"
   exit 1 # Exit script after printing help
}

while getopts "F:R:" opt
do
   case "$opt" in
      F ) parameterF="$OPTARG" ;;
      R ) parameterR="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$parameterF" ] || [ -z "$parameterR" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

src=$PWD

mkdir $src/raw_fastqc
mkdir $src/raw_fastqc/fastqc_F

for f in $src/raw_data/*$parameterF;
do fastqc $f -o $src/raw_fastqc/fastqc_F
done

multiqc $src/raw_fastqc/fastqc_F -o $src/raw_fastqc/fastqc_F/multiqc


mkdir $src/raw_fastqc/fastqc_R

for f in $src/raw_data/*$parameterR;
do fastqc $f -o $src/raw_fastqc/fastqc_R
done

multiqc $src/raw_fastqc/fastqc_R -o $src/raw_fastqc/fastqc_R/multiqc

mkdir $src/quality_reports
mv $src/raw_fastqc/fastqc_F/multiqc/multiqc_report.html $src/quality_reports/Raw_data_R1_multiqc_report.html
mv $src/raw_fastqc/fastqc_R/multiqc/multiqc_report.html	$src/quality_reports/Raw_data_R2_multiqc_report.html

