#!/bin/bash

set -eu

proj=u2023008
mail=kristina.benevides@umu.se
in=/mnt/ada/projects/aspseq/mschmid/poplar-CRISPR-WGS/RNA-seq/alternative_splicing/MEME
infile=Line3_IR_donor_sequences.fa
out=$in/pictogram_logos/Line3_donor
sif=$(realpath ../singularity/memesuite_5.5.2.sif)

if [ ! -d $out ]; then
	mkdir -p $out
fi

export APPTAINER_BINDPATH="/mnt:/mnt"

sbatch --mail-user=$mail -A $proj -o ../logs/meme_logo.out -e ../logs/meme_logo.err \
runMEMEMakeLogo.sh $sif $in/$infile $out
