#!/bin/bash -l
#SBATCH -p core
#SBATCH -t 2-00:00:00
#SBATCH -n 20
#SBATCH -A u2023008

set -eu
export SINGULARITY_BINDPATH="/mnt:/mnt"

bcftools=/mnt/picea/storage/singularity/kogia/bcftools_1.16.sif
ref=/mnt/ada/projects/aspseq/mschmid/poplar-CRISPR-WGS/WGS/reference/primary-plus-alternative-haplotypes.fasta

cd /mnt/ada/projects/aspseq/mschmid/poplar-CRISPR-WGS/WGS/bcftoolstest/
$bcftools mpileup -Ou -f $ref -b bamlist.txt --min-MQ 10 --min-BQ 20 --max-depth 1000 | $bcftools call -mv -Ob -o calls.bcf
