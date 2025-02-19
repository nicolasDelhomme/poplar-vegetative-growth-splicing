#!/bin/bash -l

set -eu

# define variables
proj=u2023008
email=kristina.benevides@umu.se
in=/mnt/ada/projects/aspseq/mschmid/poplar-CRISPR-WGS/RNA-seq/preprocessed/STAR
out=/mnt/ada/projects/aspseq/mschmid/poplar-CRISPR-WGS/RNA-seq/preprocessed/STAR/merged
reference=$(realpath ../reference)
samtools_sif=$(realpath ../singularity/kogia/samtools_1.16.sif)

if [ ! -d $out ]; then
  mkdir -p $out
fi

# env
export SINGULARITY_BINDPATH="/mnt:/mnt"

prefixes=("L03" "L26" "T89")

for prefix in "${prefixes[@]}"; do
    input_files=("$in/${prefix}_*".bam)  # Collect files with the same prefix
    sbatch -A $proj --mail-user $email -J "${prefix}_BAM_merge" \
         runMergeBams.sh $samtools_sif $input_files $out
done

