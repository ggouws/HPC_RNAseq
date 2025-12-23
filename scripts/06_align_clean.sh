#!/bin/bash

#SBATCH --job-name=06_align_clean
#SBATCH -o 06_align_clean_o%j
#SBATCH -e 06_align_clean_e%j
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4GB
#SBATCH --time=48:00:00

source ~/.bash_profile
conda activate rnaseq

src=$PWD

# make new directories for flagstat reports and cleaned bam files
mkdir $src/aligned_clean

# run flagstat to check mapping efficiency 
# use samtools to include reads mapped in pairs (-f 2) and remove all unmapped or single mapped reads (-F 12) from the bam file
# cleaned bam files are indexed (not needed for subsequent steps but handy if data are to be downloaded and viewed in IGV, for example).
# rerun flagstat to check read numbers in clean bam files

for f in $src/aligned/*.bam;
do
FBASE=$(basename $f);
BASE=${FBASE%.bam};
samtools view -q 40 -f 2 -F 12 -b $src/aligned/${BASE}.bam > $src/aligned_clean/${BASE}.clean.bam;
samtools index $src/aligned_clean/${BASE}.clean.bam;
echo ${BASE} >> $src/quality_reports/final_mapping_quality;
echo "" >> $src/quality_reports/final_mapping_quality;
samtools flagstat $src/aligned_clean/${BASE}.clean.bam >> $src/quality_reports/final_mapping_quality;
echo "" >> $src/quality_reports/final_mapping_quality;
echo "" >> $src/quality_reports/final_mapping_quality;
done


