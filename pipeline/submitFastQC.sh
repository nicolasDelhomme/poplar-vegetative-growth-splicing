#!/bin/bash -l

set -eu

# define variables
proj=u2023008
email=kristina.benevides@umu.se
in=$(realpath ../data/WGS)
out=$(realpath ../data/fastqc)
fastqc_sif=$(realpath ../singularity/kogia/fastqc_0.11.9.sif)

if [ ! -d $out ]; then
  mkdir -p $out
fi

# env
export SINGULARITY_BINDPATH="/mnt:/mnt"

# run
for f in $(find $in -name "*.fq.gz"); do

sbatch -A $proj -p core --mail-user $email -J $(basename ${f/_1.fq.gz/}) \
	runFastQC.sh $fastqc_sif $out $f
done
