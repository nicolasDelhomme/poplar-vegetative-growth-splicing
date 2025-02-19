#!/bin/bash -l
#SBATCH --mail-type=all
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00

# safeguard
set -eu

# functions
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# usage
USAGETXT=\
"
    Usage: runBwamem.sh [options] <singularity bwa container> <singularity samtools container> <Index> <Output dir> <Forward read> <Reverse read>

    Optional arguments:
    -t <INT>     ----- Number of mapping and sorting threads
    -I           ----- Index the fasta file if no index is present
    -O <string>  ----- optional arguments to pass to bwa mem (must be quoted)
    -s           ----- single end alignment
"

# VARS
OPTS=""
SINGLE=0
DOINDEX=0
THREADS=1
REV=

# options
while getopts t:O:Is opt
do
  case "$opt" in
	t) THREADS=$OPTARG;;
	I) DOINDEX=1;;
	O) OPTS=$OPTARG;;
	s) SINGLE=1;;
	\?) usage;;
  esac
done

shift `expr $OPTIND - 1`

# tests
ARGSN=6
[[ $SINGLE -eq 1 ]] && ARGSN=5

[[ $# -ne $ARGSN ]] && abort "The script expects $ARGSN arguments"

[[ ! -f $1 ]] && abort "The first argument needs to be a file, the bwa singularity container"

[[ ! -f $2 ]] && abort "The second argument needs to be a file, the samtools singularity container"

[[ ! -f $3 ]] && abort "The third argument needs to be a file, the bwa index"

[[ ! -d $4 ]] && abort "The fourth argument needs to be a directory"

[[ ! -f $5 ]] && abort "The fifth argument needs to be a file, the forward sequencing reads file"

[[ $ARGSN -eq 6 ]] && [[ ! -f $6 ]] && abort "The sixth argument needs to be a file, the reverse sequencing reads file"

[[ ${SINGULARITY_BINDPATH:-1} -eq 1 ]] && abort "This function relies on singularity, set the SINGULARITY_BINDPATH environment variable"

# index if set and needed
if [ $DOINDEX -eq 1 ]; then
    NEEDS_INDEX=0
    for suffix in .ann .pac .bwt .sa .amb .fai; do
	    if [ ! -f ${INDEX}$suffix ]; then
	      NEEDS_INDEX=1
	    fi
    done
    if [ NEEDS_INDEX ]; then
	    singularity exec $1 bwa index $3
    fi
fi

# run
[[ $ARGSN -eq 6 ]] && REV=$6
fnam=$4/$(basename ${5/_1.fq.gz/}).sorted.bam


singularity exec $1 bwa mem $OPTS -t $THREADS $3 $5 $REV | \
singularity exec $2 samtools view -bT $3 - | \
singularity exec $2 samtools sort -@ $THREADS -o $fnam

singularity exec $2 samtools index $fnam
