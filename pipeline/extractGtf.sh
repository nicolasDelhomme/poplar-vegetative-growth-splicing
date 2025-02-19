#!/bin/bash

# GFF file path
gff_file="/mnt/picea/storage/reference/Populus-tremula/v2.2/gff/Potra02_genes.gff"

# Gene IDs file path
gene_ids_file="../doc/SM_LSM_gene_ids.txt"

# Output file
output_file="../results/SM_LSM_genes.gff"

# Loop through gene IDs
while IFS= read -r gene_id; do
    echo "Processing gene ID: $gene_id"
    awk -v gene_id="$gene_id" '($3 == "gene" && $9 ~ gene_id)' "$gff_file" >> "$output_file"
done < "$gene_ids_file"

