# The splicing genes SmEa and SmEb regulate plant development during vegetative growth in poplar

__Goretti D.1,2&, Collani S.1,2, Marcon A.3, Nilsson O.3, Schmid M.1,4&__

1. Umeå Plant Science Centre, Department of Plant Physiology, Umeå University, SE-90187 Umeå, Sweden
2. DISSTE, University of Eastern Piedmont, Vercelli, 13100, Italy
3. Umeå Plant Science Centre, Department of Forest Genetics and Plant Physiology, Swedish University of Agricultural Sciences, SE-90183 Umeå, Sweden
4. Department of Plant Biology, Linnean Center for Plant Biology, Swedish University of Agricultural Sciences, S-75007 Uppsala, Sweden

& corresponding author

This repository DOI: [![DOI](https://zenodo.org/badge/935348795.svg)](https://doi.org/10.5281/zenodo.14892044)

The original repository DOI: [![DOI](https://zenodo.org/badge/706585581.svg)](https://doi.org/10.5281/zenodo.14892921)

## Abstract

Spliceosomes, the molecular machinery that catalyzes RNA splicing, contain at their core heptameric rings of Sm (or LSm) proteins. Loss of Sm gene function can have detrimental consequences on plant growth and development, as exemplified by the strong temperature-sensitive phenotype of the sme1/pcp mutant in the model plant species Arabidopsis thaliana. In this work, we characterized the SmE genes from poplar. We found that poplar contains two paralogs SmE genes, PtSmEa and PtSmEb, that encode for identical proteins. Mutagenesis using CRISPR/Cas9 revealed both overlapping and specific functions of the two poplar SmE genes, which was further substantiated by transcriptome and alternative splicing analyses using RNA-seq. Our study provides a first glimpse at the function of Sm genes in a woody species and suggests a key role of PtSmEs in regulating vegetative growth in trees.

## Setup

```bash
# data
mkdir data
ln -s /mnt/ada/projects/aspseq/mschmid/poplar-CRISPR-WGS/RNA-seq/raw data/RNASeq
ln -s /mnt/ada/projects/aspseq/mschmid/poplar-CRISPR-WGS/raw data/WGS

# tools
ln -s /mnt/picea/storage/singularity .

# reference
mkdir reference
ln -s /mnt/picea/storage/reference/rRNA/sortmerna/v4.3.4 reference/sortmerna
ln -s /mnt/picea/storage/reference/Populus-tremula/v2.2/indices/salmon1.6.0/ reference/salmon
ln -s /mnt/picea/storage/reference/Populus-tremula/v2.2/annotation reference/annotation
ln -s /mnt/picea/storage/reference/Illumina/adapters/TruSeq3-PE-2.fa reference/trimmomatic

# results
ln -s /mnt/ada/projects/aspseq/mschmid/poplar-CRISPR-WGS/RNA-seq/preprocessed analysis
```

## Description

There are two CRISPR lines with different targets:

* For the Line 26, the target locus is `Potri.006G174000/Potra001877g14982` (Potra2n6c13821)
* For the Line 03 it is `Potri.018G096200/Potra001273g10998` ( 	Potra2n6c13821)

## Note

Because of ownership issues, the original repository could not be initially made public, hence the two repositories with different DOIs but similar content.
