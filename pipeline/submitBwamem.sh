#!/bin/bash -l

set -eu

# define variables
proj=u2023008
email=kristina.benevides@umu.se
in=$(realpath ../data/WGS/raw/files)
out=/mnt/ada/projects/aspseq/mschmid/poplar-CRISPR-WGS/bwa/T89
bwa_sif=$(realpath ../singularity/kogia/bwa_0.7.17.sif)
samtools_sif=$(realpath ../singularity/kogia/samtools_1.16.sif)
index=/mnt/picea/storage/reference/Populus-tremula_X_Populus-tremuloides/v2.0/indices/bwa/primary-plus-alternative-haplotypes.fasta.gz

if [ ! -d $out ]; then
  mkdir -p $out
fi

# env
export SINGULARITY_BINDPATH="/mnt:/mnt"

# run
for f in $(find $in -name "*_1.fq.gz"); do

ID=$(basename ${f/_1.fq.gz/} | cut -f1,2,3 -d"_")
PU=$(basename ${f/_1.fq.gz/} | cut -f4,5 -d"_" | sed 's/L//')
SM=$(basename ${f/_1.fq.gz/} | cut -f1,2 -d"_")

sbatch -A $proj --mail-user $email -J $(basename ${f/_1.fq.gz/}) -o ../logs/'bwa-mem_'$(basename ${f/_1.fq.gz/}).out \
       -e ../logs/'bwa-mem_'$(basename ${f/_1.fq.gz/}).err \
       runBwamem.sh -O "-R "'@RG\tID:'${ID}'\tLB:WGS-'${SM}'\tPL:Illumina\tPU:'${PU}'\tSM:'${SM}"" \
       -t 12 $bwa_sif $samtools_sif $index $out $f ${f/_1.fq.gz/_2.fq.gz}
done
