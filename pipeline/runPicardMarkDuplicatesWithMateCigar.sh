#!/bin/bash -l

#SBATCH -t 4:00:00
#SBATCH -p core
#SBATCH -n 1
#SBATCH --mem=16G
#SBATCH --error=../logs/picard.%J.err
#SBATCH --output=../logs/picard.%J.out

# safeguard
set -eu

# functions
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# usage
USAGETXT=\
"
  Usage: $0 [options] <singularity picard container> <singularity samtools container> <BAM file> <output directory>
"
# vars
THREADS=1
JavaMem=16G
MIN=1000

# options
while getopts j:m: option
do
    case "$option" in
	    j) JavaMem=$OPTARG;;
      m) MIN=$OPTARG;;
	    \?) usage;;
    esac
done
shift `expr $OPTIND - 1`

[[ $# -ne 4 ]] && abort "This script expects three arguments"

[[ ! -f $1 ]] && "The first argument needs to be a file, the picard singularity container"

[[ ! -f $2 ]] && "The second argument needs to be a file, the picard singularity container"

[[ ! -f $3 ]] && "The third argument needs to be a bam file"

[[ ! -d $4 ]] && "The fourth argument needs to be an existing directory file"

[ -z ${SINGULARITY_BINDPATH:-} ] && abort "This function relies on singularity, set the SINGULARITY_BINDPATH environment variable"

# output
fnam=$4/$(basename ${3/.bam/})

# Run MarkDuplicates
singularity exec $1 java -Xmx${JavaMem} -XX:ParallelGCThreads=$THREADS \
    -jar /usr/picard/picard.jar MarkDuplicatesWithMateCigar \
    ASSUME_SORTED=true \
    INPUT=$3 \
    OUTPUT=${fnam}_mkdup.bam \
    METRICS_FILE=${fnam}_mkdup.metrics \
    VALIDATION_STRINGENCY=LENIENT \
    MINIMUM_DISTANCE=$MIN

# Index
singularity exec $2 samtools index ${fnam}_mkdup.bam
