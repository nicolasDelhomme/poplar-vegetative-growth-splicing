#!/bin/bash -l
#SBATCH -p core
#SBATCH -t 1-00:00:00
#SBATCH -n 2
#SBATCH -A u2023008

# extract region arounf GOI for re-alignment with STAR on T89
# GOIs
# Name	V2
# SNRPG	Potra2n6c13472
# SNRPB/N	Potra2n11c22477
# SNRPD1	Potra2n5c10994
# SNRPD2	Potra2n14c27408
# SNRPD3	Potra2n2c6364
# SNRPF	Potra2n8c17320
# SNRPE	Potra2n6c13821
# LSM1A-B	Potra2n1c1466
# LSM2	Potra2n4c10338
# LSM3A-B	Potra2n5c11140
# LSM4	Potra2n2c4537
# LSM6A-B	Potra2n8c17320
# LSM7	Potra2n1c2324
# LSM8	Potra2n1c2419
# grep -f goiID.txt /mnt/picea/storage/reference/Populus-tremula/v2.2/gff/Potra02_genes.gff | awk '{if($3 == "gene") print $1 "\t" $3 "\t" $4-1000 "\t" $5+1000 "\t" $9;}' | sort
# chr11	gene	1122223	1128046	ID=Potra2n11c22477;Name=Potra2n11c22477
# chr14	gene	8709596	8715019	ID=Potra2n14c27408;Name=Potra2n14c27408
# chr1	gene	14807932	14813592	ID=Potra2n1c1466;Name=Potra2n1c1466
# chr1	gene	27657206	27662706	ID=Potra2n1c2324;Name=Potra2n1c2324
# chr1	gene	29010878	29019018	ID=Potra2n1c2419;Name=Potra2n1c2419
# chr2	gene	22750009	22754243	ID=Potra2n2c6364;Name=Potra2n2c6364
# chr2	gene	7557288	7561879	ID=Potra2n2c4537;Name=Potra2n2c4537
# chr4	gene	19352968	19358355	ID=Potra2n4c10338;Name=Potra2n4c10338
# chr5	gene	3546983	3552147	ID=Potra2n5c10994;Name=Potra2n5c10994
# chr5	gene	4876248	4881093	ID=Potra2n5c11140;Name=Potra2n5c11140
# chr6	gene	5162267	5166840	ID=Potra2n6c13472;Name=Potra2n6c13472
# chr6	gene	8711854	8718644	ID=Potra2n6c13821;Name=Potra2n6c13821
# chr8	gene	4637051	4642993	ID=Potra2n8c17320;Name=Potra2n8c17320
# grep 'Potra2n18c32411' /mnt/picea/storage/reference/Populus-tremula/v2.2/gff/Potra02_genes.gff | awk '{if($3 == "gene") print $1 "\t" $3 "\t" $4-1000 "\t" $5+1000 "\t" $9;}' | sort
# chr18   gene    3716950 3723446 ID=Potra2n18c32411;Name=Potra2n18c32411

set -eu
export SINGULARITY_BINDPATH="/mnt:/mnt"

in=$(realpath ../analysis/STAR)
out=$(realpath ../analysis/T89_STAR/GOIreads)
samtools=$(realpath ../singularity/kogia/samtools_1.16.sif)

mkdir -p $out

# run
for f in $(find $in -name "*_Aligned.sortedByCoord.out.bam"); do
$samtools index $f
$samtools view -hb $f \
"chr1:14807932-14813592" \
"chr1:27657206-27662706" \
"chr1:29010878-29019018" \
"chr2:22750009-22754243" \
"chr2:7557288-7561879" \
"chr4:19352968-19358355" \
"chr5:3546983-3552147" \
"chr5:4876248-4881093" \
"chr6:5162267-5166840" \
"chr6:8711854-8718644" \
"chr8:4637051-4642993" \
"chr11:1122223-1128046" \
"chr14:8709596-8715019" \
"chr18:3716950-3723446" \
| $samtools sort -n -o $out/$(basename ${f/_Aligned.sortedByCoord.out.bam/})_goireads.bam -
$samtools fastq -1 $out/$(basename ${f/_Aligned.sortedByCoord.out.bam/})_goireads_1.fq.gz -2 $out/$(basename ${f/_Aligned.sortedByCoord.out.bam/})_goireads_2.fq.gz $out/$(basename ${f/_Aligned.sortedByCoord.out.bam/})_goireads.bam
done
