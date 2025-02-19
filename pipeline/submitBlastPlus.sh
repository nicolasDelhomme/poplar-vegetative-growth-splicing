#!/bin/bash -l

## error and verbose
set -ex

## default args
account=u2023008
mail="kristina.benevides@umu.se"
sif=$(realpath ../singularity/ncbi-blast_2.13.0.sif)
in=$(realpath ../results/SM_LSM_genes.fasta)
out=/mnt/ada/projects/aspseq/mschmid/poplar-CRISPR-WGS/WGS/blast
inx=/mnt/ada/projects/aspseq/sjansson/T89_BLAST/blast_db/primary-plus-alternative-haplotypes.fasta.gz
#inx=/mnt/picea/storage/reference/Populus-tremula/v2.2/indices/BLAST+/Potra02_genome.fasta
## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

## env
export SINGULARITY_BINDPATH="/mnt:/mnt"

## prepare
sbatch -A $account -p core --mail-user $mail -e $out/BlastPlus.err -o $out/BlastPlus.out \
runBlastPlus.sh $sif blastn $in $inx $out

