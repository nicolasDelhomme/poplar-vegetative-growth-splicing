#!/bin/bash -l

#safeguard

set -eu

# define variables
proj=u2023008
email=kristina.benevides@umu.se
in=/mnt/ada/projects/aspseq/mschmid/poplar-CRISPR-WGS/WGS/picard/T89
out=/mnt/ada/projects/aspseq/mschmid/poplar-CRISPR-WGS/WGS/gatk/T89
ref=$(realpath ../reference/v2.0/fasta/primary-plus-alternative-haplotypes.fasta.gz)
gatk_sif=$(realpath ../singularity/gatk_4.2.6.1.sif)
dbsnp_dir=/mnt/ada/projects/aspseq/mschmid/poplar-CRISPR-WGS/WGS/gatk/known_sites

if [ ! -d $out ]; then
  mkdir -p $out
fi

# env
export SINGULARITY_BINDPATH="/mnt:/mnt"

# run
for BAM in $(find $in -name "*.sorted_mkdup.bam"); do
bname=$(basename ${BAM/.sorted_mkdup.bam/})
sbatch -A $proj --mail-user $email -J 'gatk.'$bname \
        runGATK.sh $gatk_sif $BAM $ref $out $dbsnp_dir/$bname'.sort.dup.filtered.hc.vcf.gz'
done
