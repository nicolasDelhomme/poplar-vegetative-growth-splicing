#!/bin/bash
#SBATCH -n 8
#SBATCH -t 2-00:00:00

# source helpers
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

USAGETXT=\
"
Usage: $(basename $0) <singularity container> <input FASTA file> <output directory>
"

## arguments
[[ $# -ne 3 ]] && abort "This script takes three arguments"

[[ ! -f $1 ]] && abort "The first argument needs to be an existing singularity container file."

## enforce apptainer
[[ ${APPTAINER_BINDPATH:-1} -eq 1 ]] && abort "This function relies on apptainer, set the APPTAINER_BINDPATH environment variable"

[[ ! -f $2 ]] && abort "The second argument needs to be an existing FASTA file to search for motifs."

[[ ! -d $3 ]] && abort "The third argument needs to be an existing output directory."

# Run MEME Suite
apptainer exec $1 meme \
    $2 \
    -oc $3 \
    -dna \
    -nmotifs 1 \
    -minw 6 \
    -maxw 50 \
    -mod anr
