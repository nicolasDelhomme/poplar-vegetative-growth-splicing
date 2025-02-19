#!/bin/bash -l

# be safe (e stops on error; u checks for undefined variables)
set -eu


# variables
start=2
end=7
proj=u2023008
mail=aman.zare@umu.se
singularity=$(realpath ../singularity/kogia)
reference=$(realpath ../reference)
sortmerna_fasta=$reference/sortmerna/fasta/smr_v4.3_fast_db.fasta
sortmerna_inx=$reference/sortmerna/indices/smr_v4.3_fast_db
trimmomatic_adapter=$reference/trimmomatic/TruSeq3-PE-2.fa
salmon_index=$reference/salmon
in=$(realpath ../data/RNASeq)
out=$(realpath ../analysis/preprocessed)

# create the out dir if it does not exist
[[ ! -d $out ]] && mkdir -p $out

#if [ ! -d $out ]; then
#  mkdir -p $out
#fi

# singularity
export SINGULARITY_BINDPATH="/mnt:/mnt"

# loop for every file in the in directory, call the preprocessing
for f in $(find $in -name "*_1.fq.gz"); do
  bash ../UPSCb-common/pipeline/runRNASeqPreprocessing.sh \
  -s $start -e $end \
  -x $sortmerna_inx -X $sortmerna_fasta \
  -A $trimmomatic_adapter \
  -S $salmon_index \
  $proj $mail $singularity $f ${f/_1.fq/_2.fq} $out 
done
