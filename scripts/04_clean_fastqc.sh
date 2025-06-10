#!/bin/bash

#SBATCH --job-name=04_clean_fastqc
#SBATCH -o 04_clean_fastqc_o%j
#SBATCH -e 04_clean_fastqc_e%j
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH -A molecolb
#SBATCH -p molecolb
#SBATCH --mem-per-cpu=12GB
#SBATCH --time=48:00:00

source ~/.bash_profile
conda activate multiqc

src=$PWD

mkdir $src/clean_fastqc
mkdir $src/clean_fastqc/fastqc_F

for f in $src/trimmed/*_trimmed_paired_R1.fastq.gz;
do fastqc $f -o $src/clean_fastqc/fastqc_F
done

multiqc $src/clean_fastqc/fastqc_F -o $src/clean_fastqc/fastqc_F/multiqc


mkdir $src/clean_fastqc/fastqc_R

for f in $src/trimmed/*_trimmed_paired_R2.fastq.gz;
do fastqc $f -o $src/clean_fastqc/fastqc_R
done

multiqc $src/clean_fastqc/fastqc_R -o $src/clean_fastqc/fastqc_R/multiqc

mv $src/clean_fastqc/fastqc_F/multiqc/multiqc_report_1.html $src/quality_reports/Clean_data_R1_multiqc_report.html
mv $src/clean_fastqc/fastqc_R/multiqc/multiqc_report_1.html $src/quality_reports/Clean_data_R2_multiqc_report.html

