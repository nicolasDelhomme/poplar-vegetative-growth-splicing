#!/bin/bash -l
#SBATCH -p node
#SBATCH -n 20
#SBATCH -t 2-00:00:00
#SBATCH -A u2023008

set -eu
export SINGULARITY_BINDPATH="/mnt:/mnt"

in=/mnt/ada/projects/aspseq/mschmid/poplar-CRISPR-WGS/RNA-seq/preprocessed/sortmerna/fqout
out=/mnt/ada/projects/aspseq/mschmid/poplar-CRISPR-WGS/RNA-seq/preprocessed/STAR
STAR=/mnt/picea/storage/singularity/kogia/star_2.7.9a.sif
reference=/mnt/picea/storage/reference/Populus-tremula/v2.2/indices/star/Potra02_STAR_v2.7.9a

mkdir -p $out
cd $out

#T89
$STAR --genomeDir $reference \
--readFilesIn $in/T89_1_EKRN230013111-1A_HVF2VDSX5_L2_1_fwd.fq.gz,$in/T89_1_EKRN230013111-1A_HVGLKDSX5_L3_1_fwd.fq.gz  $in/T89_1_EKRN230013111-1A_HVF2VDSX5_L2_1_rev.fq.gz,$in/T89_1_EKRN230013111-1A_HVGLKDSX5_L3_1_rev.fq.gz \
--outFileNamePrefix T89_1_ \
--readFilesCommand zcat --alignIntronMax 11000 --runThreadN 20 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000

$STAR --genomeDir $reference \
--readFilesIn $in/T89_2_EKRN230013112-1A_HVGLKDSX5_L3_1_fwd.fq.gz  $in/T89_2_EKRN230013112-1A_HVGLKDSX5_L3_1_rev.fq.gz \
--outFileNamePrefix T89_2_ \
--readFilesCommand zcat --alignIntronMax 11000 --runThreadN 20 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000

$STAR --genomeDir $reference \
--readFilesIn $in/T89_3_EKRN230013113-1A_HVF2VDSX5_L2_1_fwd.fq.gz,$in/T89_3_EKRN230013113-1A_HVGLKDSX5_L3_1_fwd.fq.gz  $in/T89_3_EKRN230013113-1A_HVF2VDSX5_L2_1_rev.fq.gz,$in/T89_3_EKRN230013113-1A_HVGLKDSX5_L3_1_rev.fq.gz \
--outFileNamePrefix T89_3_ \
--readFilesCommand zcat --alignIntronMax 11000 --runThreadN 20 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000

$STAR --genomeDir $reference \
--readFilesIn $in/T89_4_EKRN230013114-1A_HVFVCDSX5_L4_1_fwd.fq.gz,$in/T89_4_EKRN230013114-1A_HVFWHDSX5_L4_1_fwd.fq.gz  $in/T89_4_EKRN230013114-1A_HVFVCDSX5_L4_1_rev.fq.gz,$in/T89_4_EKRN230013114-1A_HVFWHDSX5_L4_1_rev.fq.gz \
--outFileNamePrefix T89_4_ \
--readFilesCommand zcat --alignIntronMax 11000 --runThreadN 20 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000

#-readFilesIn PCa10_S6_L001_R1_001.fastq.gz,PCa10_S7_L001_R1_001.fastq.gz PCa10_S6_L001_R2_001.fastq.gz,PCa10_S7_L001_R2_001.fastq.g
# --readFilesIn A_R1,B_R1,C_R1  A_R2,B_R2,C_R2
#  --outSAMattrRGline:
#ID:sampleA CN:AA DS:AAA , ID:sampleBB CN:bb DS:bbbb , ID:sampleC CN:ccc DS:cccc 
#We should also be careful with the alignments. Maybe we should not allow multi-mapping. We should also run the alignments on the genome from Kristina that has both haplotypes.#