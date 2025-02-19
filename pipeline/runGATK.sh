#!/bin/bash -l

#SBATCH -t 4-00:00:00
#SBATCH -p core
#SBATCH -n 8
#SBATCH --mem 64G
#SBATCH --error=../logs/gatk.%J.err
#SBATCH --output=../logs/gatk.%J.out

# safeguard
set -eu

# helper
source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

# usage
USAGETXT=\
"
Usage: runGatk_BaseRecalibration.sh <singularity gatk container> <BAM file> <fasta ref> <output directory> <dbsnp>

Note: This script is GATK v4 compatible and V3 incompatible. More at https://software.broadinstitute.org/gatk/blog?id=7847
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

if [ ! -f $5 ]; then
    abort "Could not find dbSNP file '$5'"
fi
db=$5

[ -z ${SINGULARITY_BINDPATH:-} ] && abort "This function relies on singularity, set the SINGULARITY_BINDPATH enviroment variable"

name_out=`basename "${inbam/.bam/.table}"`
outname=`basename "${inbam/.bam/_recalibrated.bam}"`
post_out=`basename "${inbam/.bam/_recalibrated.table}"`
plots_out=`basename "${inbam/.bam/.pdf}"`
gvcf_out=`basename "${inbam/.bam/.vcf.gz}"`
vcf_out=`basename "${inbam/.bam/.snps.indels.vcf.gz}"`
filt_vcf_out="${inbam/.bam/.filtered.snps.indels.vcf.gz}"

# Run BaseRecalibrator
# Analyze patterns of covariation in the sequence dataset
singularity exec $sif gatk BaseRecalibrator -R $ref --known-sites $db -I $inbam -O $outdir/$name_out

# Apply the recalibration to your sequence data
singularity exec $sif gatk ApplyBQSR -R $ref -I $inbam --bqsr-recal-file $outdir/$name_out -O $outdir/$outname

# Do a second pass to analyze covariation remaining after recalibration
singularity exec $sif gatk BaseRecalibrator -R $ref --known-sites $db -I $outdir/$outname -O $outdir/$post_out

# Generate before/after plots
singularity exec $sif gatk AnalyzeCovariates -before $outdir/$name_out -after $outdir/$post_out -plots $outdir/$plots_out

# Run HaplotypeCaller
singularity exec $sif gatk --java-options "-Xmx40G -Djava.io.tmpdir=/mnt/picea/tmp" HaplotypeCaller \
-I $outdir/$outname -O $outdir/$vcf_out -R $ref --annotation FisherStrand --annotation QualByDepth --annotation MappingQuality \
  --annotation-group StandardAnnotation --create-output-variant-index true

# Run VariantFiltration
singularity exec $sif gatk --java-options "-Xmx40G -Djava.io.tmpdir=/mnt/picea/tmp" VariantFiltration -OVI true \
  --variant $outdir/$vcf_out \
  --output $filt_vcf_out \
  --filter-name FisherStrand --filter 'FS > 60.0' \
  --filter-name QualByDepth --filter 'QD < 2.0' \
  --filter-name MappingQuality --filter 'MQ < 40.0'
