#!/bin/bash -l
#SBATCH --mail-type=all
#SBATCH --mem=64G
#SBATCH --time=12:00:00

# safeguard
set -eu

# functions
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# usage
USAGETXT=\
"
    Usage: runBwamem.sh [options] <singularity samtools container> <BAMs> <Output Dir>

    Optional arguments:
    -O <string>  ----- optional arguments to pass to samtools (must be quoted)
"

# VARS
OPTS=""
THREADS=1

# options
while getopts O opt
do
  case "$opt" in
	O) OPTS=$OPTARG;;
	\?) usage;;
  esac
done

shift `expr $OPTIND - 1`

# tests
[[ ! -f $1 ]] && abort "The first argument needs to be the samtools singularity container"

[[ ${SINGULARITY_BINDPATH:-1} -eq 1 ]] && abort "This function relies on singularity, set the SINGULARITY_BINDPATH environment variable"

# run
fnam=$(basename ${2/.bam/.merged.bam})

singularity exec $1 samtools merge --threads 4 -o ${!#}/$fnam $2
singularity exec $1 samtools index ${!#}/$fnam
