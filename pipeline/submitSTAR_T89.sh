#!/bin/bash -l
#SBATCH -p node
#SBATCH -n 20
#SBATCH -t 2-00:00:00
#SBATCH -A u2023008

# be safe (e stops on error; u checks for undefined variables)
set -eu

# vars
in=$(realpath ../analysis/T89_STAR/GOIreads)
out=$(realpath ../analysis/T89_STAR/results)
reference=$(realpath ../analysis/T89_STAR/index)
STAR=/mnt/picea/storage/singularity/kogia/star_2.7.9a.sif

# env
export SINGULARITY_BINDPATH="/mnt:/mnt"

mkdir -p $out
cd $out

# run
for f in $(find $in -name "*_1.fq.gz"); do
$STAR --genomeDir $reference \
--readFilesIn $f ${f/_1.fq.gz/_2.fq.gz} \
--outFileNamePrefix $(basename ${f/_1.fq.gz/}) \
--readFilesCommand zcat --alignIntronMax 11000 --runThreadN 20 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000
done

# $STAR --genomeDir $reference \
# --readFilesIn $in/T89_1_2goireads_1.fq.gz $in/T89_1_2goireads_2.fq.gz \
# --outFileNamePrefix T89_1_2goireads \
# --readFilesCommand zcat --alignIntronMax 11000 --runThreadN 20 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000
# 
# $STAR --genomeDir $reference \
# --readFilesIn $in/T89_2_2goireads_1.fq.gz $in/T89_2_2goireads_2.fq.gz \
# --outFileNamePrefix T89_2_2goireads \
# --readFilesCommand zcat --alignIntronMax 11000 --runThreadN 20 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000
# 
# $STAR --genomeDir $reference \
# --readFilesIn $in/T89_3_2goireads_1.fq.gz $in/T89_3_2goireads_2.fq.gz \
# --outFileNamePrefix T89_3_2goireads \
# --readFilesCommand zcat --alignIntronMax 11000 --runThreadN 20 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000
# 
# $STAR --genomeDir $reference \
# --readFilesIn $in/L03_1_2goireads_1.fq.gz $in/L03_1_2goireads_2.fq.gz \
# --outFileNamePrefix L03_1_2goireads \
# --readFilesCommand zcat --alignIntronMax 11000 --runThreadN 20 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000

# $STAR --genomeDir $reference \
# --readFilesIn $in/L03_3_2goireads_1.fq.gz $in/L03_3_2goireads_2.fq.gz \
# --outFileNamePrefix L03_3_2goireads \
# --readFilesCommand zcat --alignIntronMax 11000 --runThreadN 20 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000

# $STAR --genomeDir $reference \
# --readFilesIn $in/L03_4_2goireads_1.fq.gz $in/L03_4_2goireads_2.fq.gz \
# --outFileNamePrefix L03_4_2goireads \
# --readFilesCommand zcat --alignIntronMax 11000 --runThreadN 20 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000

# $STAR --genomeDir $reference \
# --readFilesIn $in/L26_1_2goireads_1.fq.gz $in/L26_1_2goireads_2.fq.gz \
# --outFileNamePrefix L26_1_2goireads \
# --readFilesCommand zcat --alignIntronMax 11000 --runThreadN 20 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000
# 
# $STAR --genomeDir $reference \
# --readFilesIn $in/L26_2_2goireads_1.fq.gz $in/L26_2_2goireads_2.fq.gz \
# --outFileNamePrefix L26_2_2goireads \
# --readFilesCommand zcat --alignIntronMax 11000 --runThreadN 20 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000
# 
# $STAR --genomeDir $reference \
# --readFilesIn $in/L26_3_2goireads_1.fq.gz $in/L26_3_2goireads_2.fq.gz \
# --outFileNamePrefix L26_3_2goireads \
# --readFilesCommand zcat --alignIntronMax 11000 --runThreadN 20 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000
# 
# $STAR --genomeDir $reference \
# --readFilesIn $in/L26_4_2goireads_1.fq.gz $in/L26_4_2goireads_2.fq.gz \
# --outFileNamePrefix L26_4_2goireads \
# --readFilesCommand zcat --alignIntronMax 11000 --runThreadN 20 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000
