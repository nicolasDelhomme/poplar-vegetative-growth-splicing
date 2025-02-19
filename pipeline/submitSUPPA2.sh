#!/bin/bash

set -eu

proj=u2023008
mail=kristina.benevides@umu.se
out=/mnt/ada/projects/aspseq/mschmid/poplar-CRISPR-WGS/RNA-seq/alternative_splicing
in=$(realpath ../reference/gtf/Potra02_genes.gtf)
suppa2_sif=$(realpath ../singularity/SUPPA2.sif)

if [ ! -d $out ]; then
	mkdir -p $out
fi

export SINGULARITY_BINDPATH="/mnt:/mnt"

sbatch --mail-user=$mail -o ../logs/suppa2.out -e ../logs/suppa2.err -A $proj runSUPPA2GenerateEvents.sh $suppa2_sif $in $out
