#!/bin/bash

set -eu

proj=u2023008
mail=kristina.benevides@umu.se
in=/mnt/ada/projects/aspseq/mschmid/poplar-CRISPR-WGS/WGS
out=$in/multiqc
multiqc_sif=$(realpath ../singularity/multiqc_1.18.sif)

if [ ! -d $out ]; then
	mkdir -p $out
fi

export SINGULARITY_BINDPATH="/mnt:/mnt"

sbatch --mail-user=$mail -o ../logs/multiqc.out -e ../logs/multiqc.err -A $proj runMultiQC.sh $multiqc_sif  $in $out
