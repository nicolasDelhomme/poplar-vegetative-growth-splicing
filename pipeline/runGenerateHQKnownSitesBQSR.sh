#!/bin/bash -l

#SBATCH -t 1-00:00:00
#SBATCH -p core
#SBATCH -n 6
#SBATCH --mem 36G
#SBATCH --error=../logs/generate_known_sites_BQSR.%J.err
#SBATCH --output=../logs/generate_known_sites_BQSR.%J.out

# safeguard
set -eu

# helper
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# usage
USAGETXT=\
"
Usage: runGenerateHQKnownSitesBQSR.sh <singularity gatk container> <BAM file> <fasta ref> <outdir>

Note: This script is GATK v4 compatible and V3 incompatible. More at https://software.broadin>
"

# Tests
if [ ! -f $1 ]; then
     abort "Could not find singularity container '$1'"
fi
sif=$1

if [ ! -f $2 ]; then
     abort "Could not find BAM file '$2'"
fi
inbam=$2

if [ ! -f $3 ]; then
     abort "Could not find FASTA file '$3'"
fi
ref=$3

if [ ! -d $4 ]; then
     abort "Could not find directory '$4'"
fi
outdir=$4

[ -z ${SINGULARITY_BINDPATH:-} ] && abort "This function relies on singularity, set the SINGULARITY_BINDPATH enviroment"

name_out_hc=`basename "${inbam/.sorted_mkdup.bam/.sort.dup.raw.hc.vcf.gz}"`
name_out_filt=`basename "${inbam/.sorted_mkdup.bam/.sort.dup.filtered.hc.vcf.gz}"`

apptainer exec $sif gatk HaplotypeCaller -OVI true \
  --annotation FisherStrand -A QualByDepth -A MappingQuality -G StandardAnnotation \
  --input $inbam --output $outdir/$name_out_hc \
  --reference $ref

apptainer exec $sif gatk VariantFiltration -OVI true \
  --variant $outdir/$name_out_hc \
  --output $outdir/$name_out_filt \
  --filter-name FisherStrand --filter 'FS > 50.0' \
  --filter-name QualByDepth --filter 'QD < 4.0' \
  --filter-name MappingQuality --filter 'MQ < 50.0'
