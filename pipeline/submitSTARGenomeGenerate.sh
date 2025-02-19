#!/bin/bash
#SBATCH -p core
#SBATCH -t 12:00:00
#SBATCH -n 4
#SBATCH -A u2023008

set -eu
export SINGULARITY_BINDPATH="/mnt:/mnt"

proj=u2023008
outdir=/mnt/ada/projects/aspseq/mschmid/poplar-CRISPR-WGS/RNA-seq/preprocessed/T89_STAR/
star=/mnt/picea/storage/singularity/kogia/star_2.7.10a.sif

mkdir -p $outdir/index
zcat /mnt/picea/storage/reference/Populus-tremula_X_Populus-tremuloides/v2.0/fasta/primary-plus-alternative-haplotypes.fasta.gz > $outdir/genome.fa
sbatch -A $proj ./runSTARGenomeGenerate.sh -n $star $outdir/index $outdir/genome.fa
