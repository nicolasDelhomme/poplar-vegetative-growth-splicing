#!/bin/bash -l
#SBATCH -p core
#SBATCH -t 2-00:00:00
#SBATCH -n 20
#SBATCH -A u2023008

set -eu
export SINGULARITY_BINDPATH="/mnt:/mnt"

ref=/mnt/ada/projects/aspseq/mschmid/poplar-CRISPR-WGS/WGS/reference/primary-plus-alternative-haplotypes.fasta

cd /mnt/ada/projects/aspseq/mschmid/poplar-CRISPR-WGS/WGS/angsd_test/test3
~/apptainerbuilds/angsd.sif -GL 2 -doGlf 2 -b bamlist.txt -doBcf 1 -doMajorMinor 4 -ref $ref -SNP_pval 1e-6 -doMaf 2 -nThreads 20
