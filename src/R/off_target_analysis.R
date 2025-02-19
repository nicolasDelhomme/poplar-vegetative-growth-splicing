#' ---
#' title: "Poplar CRISPR lines WGS off-target analysis"
#' author: "Kristina Benevides & Fai Theerarat Kochakarn"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    fig_width: 9
#'    fig_height: 6
#'    toc: true
#'    number_sections: true
#'    toc_depth: 4
#'    toc_float:
#'      collapsed: TRUE
#'      smooth_scroll: TRUE
#'    code_folding: hide
#'    theme: "flatly"
#'    highlight: pygments
#'    includes:
#'      before_body: header.html
#'      after_body: footer.html
#'    css: style.css
#' ---
#' 
#' <hr />
#' &nbsp;
#' 

#' # Setup

#' * Libraries
suppressPackageStartupMessages({
  library(data.table)
  library(splitstackshape)
  library(tidyverse)
  library(VariantAnnotation)
  library(RColorBrewer)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  library(R.utils)
  library(Biostrings)
})


#' * Helper functions
#' 
#' Basic filtering
QC_filtering <- function(data_depth, data_qual) {
  x <- as.data.frame(table(data_depth$SUM_DEPTH))
  lower <- 0.75 * median(data_depth$SUM_DEPTH)
  upper <- median(data_depth$SUM_DEPTH) + 1 * sd(data_depth$SUM_DEPTH)
  xupper <- ceiling(upper/100) * 100
  g1 <- ggplot(x, aes(x = as.numeric(as.character(Var1)), y = Freq)) + geom_line() + xlab("Depth") +
    ylab("bp") + xlim(0, xupper) + geom_vline(xintercept = lower, color = "red",
                                              linewidth = 1.3) + geom_vline(xintercept = upper, color = "red", linewidth = 1.3) + ggtitle("Example threshold: 0.8X median depth, median depth + 2sd") %>% suppressMessages()
  #' ## Variant quality distribution

  g2 <- ggplot(subset(data_qual, QUAL < 1000), aes(x = QUAL)) + geom_histogram(fill = "white",
                                                                    color = "black", bins = 50) + xlab("Quality value") + ylab("Count") + geom_vline(xintercept = 30,                                                                                                                                                     
                                                                    color = "red", linewidth = 1.3) + ggtitle("Example threshold: Q30") %>% suppressMessages()
  
  # Combine both plots into a single figure
  combined_plot <- grid.arrange(g1, g2, ncol = 2) %>% suppressMessages()
  
  print(paste("Lower DP:", lower," ","Upper DP:", upper))
  # Return the combined plot
  return(combined_plot)
}

#' Extract DP (sequencing depth) from VCF
extract_DP <- function(column) {
  dp_entries <- str_extract_all(column, "DP=\\d+")
  return(ifelse(dp_entries != "character(0)", dp_entries, ""))
}

#' Function to filter a VCF file from GATK Joint Genotyping workflow. 
#' The input is required to have the following structure: 
#' #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT sample1 sample2 sample3 ...

filterVCF <- function(vcf,file_prefix,DP_lower,DP_upper,sample_names) {
  
  #'  # Load VCF file into a data table
  vcf_file <- fread(vcf, sep="\t", skip = "#CHROM")
  
  #'  # Split the columns with multiple information fields into sepearate columns
  snps_split <- cSplit(vcf_file, 8, sep = ";", "wide")
  
  #'  # Apply the function to INFO_4 and INFO_5 columns (either of which contain the DP information) and create a new column DP
  
  snps_split <- snps_split %>%
    mutate(DP_4 = extract_DP(INFO_04),
           DP_5 = extract_DP(INFO_05)) %>% 
    unite(DP,DP_4,DP_5) %>%
    mutate(DP=str_replace_all(DP,"_","")) 
  
  snps <- as_tibble(snps_split) %>% dplyr::select(., `#CHROM`, POS, REF, ALT, QUAL, INFO_01, INFO_02, INFO_03, DP, any_of(sample_names))
  
  #'  # Rename the new columns
  snps <- snps %>% rename_with(~gsub("INFO_01", "AC", .), ends_with("_01")) %>%
    rename_with(~gsub("INFO_02", "AF", .), ends_with("_02")) %>%
    rename_with(~gsub("INFO_03", "AN", .), ends_with("_03")) %>%
    mutate(AC = str_replace(AC, "AC=", "")) %>%
    mutate(AF = str_replace(AF, "AF=", "")) %>%
    mutate(AN = str_replace(AN, "AN=", "")) %>%
    mutate(DP = str_replace(DP, "DP=", ""))
  
  snps$DP <- as.numeric(snps$DP)
  snps$AF <- as.numeric(snps$AF)
  
  snps_filtered <- snps %>% filter(DP >= DP_lower) %>% filter(DP <= DP_upper)
  snps_filtered_homo <- snps %>% filter(DP >= DP_lower) %>% filter(DP <= DP_upper) %>% filter(AF == 1.0)
  
  #'  # Write output to tsv files
  write_tsv(snps,paste0(file_prefix, "_all_SNPs.tsv"))
  write_tsv(snps_filtered,paste0(file_prefix, "_filtered_SNPs.tsv"))
  write_tsv(snps_filtered_homo,paste0(file_prefix, "_filtered_homozygous_SNPs.tsv"))
  
  return(snps_filtered_homo) 
}

#' Function to extract +/- 50bp of the SNPs/Indels
expand_sequences <- function(chrom, position, ref_genome, window_size = 50) {
  seq_start <- max(1, position - window_size)
  seq_end <- min(length(ref_genome[[chrom]]), position + window_size)
  return(cbind(seq_start,seq_end))
}

#' ## Data QC and filtering 

#' Load the T89 reference genome
ref_genome <- readDNAStringSet("reference/v2.0/fasta/primary-plus-alternative-haplotypes.fasta.gz")
ref_genome@ranges@NAMES <- sub("\\s.*", "", ref_genome@ranges@NAMES)

#' ### Line 3

#' Load VCF 
vcf <- "data/WGS/gatk/T89/gDNA_line3_EKDN230011266-1A_HVGLKDSX5_L1.sorted_mkdup.filtered.snps.indels.vcf.gz"
system(paste("srun -A u2023008 apptainer exec -B /mnt:/mnt /mnt/picea/storage/singularity/vcftools_v0.1.16.sif vcftools --gzvcf", vcf,"--site-depth 2>/dev/null"))
data_depth <- read.table("out.ldepth", header = TRUE)
system(paste("srun -A u2023008 apptainer exec -B /mnt:/mnt /mnt/picea/storage/singularity/vcftools_v0.1.16.sif vcftools --gzvcf", vcf,"--site-qual 2>/dev/null"))
data_qual <- read.table("out.lqual", header = TRUE)

#' Basic variant filtering 
QC_filtering(data_depth,data_qual)

#' Filter for high-quality SNPs/Indels
snps_line3 <- filterVCF(vcf,"results/gDNA_line3",round(lower,0),round(upper,0),c("gDNA_line3", "gDNA_line26"))

#' Extract chrom and position to a tibble
variants_line3 <- tibble(chrom=snps_line3$`#CHROM`, pos=snps_line3$POS)

#' Expand +/- 50bp around the SNPs/Indels
expanded_sequences_line3 <- variants_line3 %>%
  rowwise() %>%
  mutate(expanded_sequences = expand_sequences(chrom, pos, ref_genome))

#' View the result
print(expanded_sequences_line3)

#' Retrieve just the start and stop position as a bed file for bedtools
bed_sequences_line3 <- tibble(expanded_sequences_line3) %>% 
  mutate(start=expanded_sequences[,1], end=expanded_sequences[,2]) %>% 
  dplyr::select(chrom, start, end)

#' Write to file
write_tsv(bed_sequences_line3, "results/snps_indels_line3.bed")

#' ### Line 26
#' Load VCF 
vcf <- "data/WGS/gatk/T89/gDNA_line26_EKDN230011267-1A_HVCGCDSX5_HYCW3DSX5_L3_L4_merged.sorted_mkdup.filtered.snps.indels.vcf.gz"
system(paste("srun -A u2023008 apptainer exec -B /mnt:/mnt /mnt/picea/storage/singularity/vcftools_v0.1.16.sif vcftools --gzvcf", vcf,"--site-depth 2>/dev/null"))
data_depth <- read.table("out.ldepth", header = TRUE)
system(paste("srun -A u2023008 apptainer exec -B /mnt:/mnt /mnt/picea/storage/singularity/vcftools_v0.1.16.sif vcftools --gzvcf", vcf,"--site-qual 2>/dev/null"))
data_qual <- read.table("out.lqual", header = TRUE)

#' Basic variant filtering 
QC_filtering(data_depth,data_qual)

#' Filter for high-quality SNPs/Indels
snps_line26 <- filterVCF(vcf,"results/gDNA_line26",round(lower,0),round(upper,0),c("gDNA_line3", "gDNA_line26"))

#' Extract chrom and position to a tibble
variants_line26 <- tibble(chrom=snps_line26$`#CHROM`, pos=snps_line26$POS)

#' Expand +/- 50bp around the SNPs/Indels
expanded_sequences_line26 <- variants_line26 %>%
  rowwise() %>%
  mutate(expanded_sequences = expand_sequences(chrom, pos, ref_genome))

#' View the result
print(expanded_sequences_line26)

#' Retrieve just the start and stop position as a bed file for bedtools
bed_sequences_line26 <- tibble(expanded_sequences_line26) %>% 
  mutate(start=expanded_sequences[,1], end=expanded_sequences[,2]) %>% 
  dplyr::select(chrom, start, end)

#' Write to file
write_tsv(bed_sequences_line26, "results/snps_indels_line26.bed")

#' ### Using Potra v2.2 as reference

ref_genome <- readDNAStringSet("reference/fasta/Potra02_genome_hardmasked.fasta.gz")

#' #### Line 3

#' Load VCF 
vcf_potra <- "data/WGS/gatk/gDNA_line3_EKDN230011266-1A_HVGLKDSX5_L1.sorted_mkdup.filtered.snps.indels.vcf.gz"
system(paste("srun -A u2023008 apptainer exec -B /mnt:/mnt /mnt/picea/storage/singularity/vcftools_v0.1.16.sif vcftools --gzvcf", vcf,"--site-depth 2>/dev/null"))
data_potra <- read.table("out.ldepth", header = TRUE)

#' Basic variant filtering 
QC_filtering(data_potra)
#' Filter for high-quality SNPs/Indels
snps_line3_potra <- filterVCF(vcf_potra,"results/gDNA_line3",round(lower,0),round(upper,0),c("gDNA_line3", "gDNA_line26"))

#' Extract chrom and position to a tibble
variants_line3_potra <- tibble(chrom=snps_line3_potra$`#CHROM`, pos=snps_line3_potra$POS)

#' Expand +/- 50bp around the SNPs/Indels
expanded_sequences_line3_potra <- variants_line3_potra %>%
  rowwise() %>%
  mutate(expanded_sequences = expand_sequences(chrom, pos, ref_genome))

#' View the result
print(expanded_sequences_line3_potra)

#' Retrieve just the start and stop position as a bed file for bedtools
bed_sequences_line3_potra <- tibble(expanded_sequences_line3_potra) %>% 
  mutate(start=expanded_sequences[,1], end=expanded_sequences[,2]) %>% 
  dplyr::select(chrom, start, end)

#' Write to file
write_tsv(bed_sequences_line3_potra, "results/snps_indels_line3_Potra02.bed")

#' #### Line 26

#' Load VCF 
vcf_potra <- "data/WGS/gatk/gDNA_line3_EKDN230011266-1A_HVGLKDSX5_L1.sorted_mkdup.filtered.snps.indels.vcf.gz"
system(paste("srun -A u2023008 apptainer exec -B /mnt:/mnt /mnt/picea/storage/singularity/vcftools_v0.1.16.sif vcftools --gzvcf", vcf,"--site-depth 2>/dev/null"))
data_potra_depth <- read.table("out.ldepth", header = TRUE)

system(paste("srun -A u2023008 apptainer exec -B /mnt:/mnt /mnt/picea/storage/singularity/vcftools_v0.1.16.sif vcftools --gzvcf", vcf,"--site-qual 2>/dev/null"))
data_potra_qual <- read.table("out.lqual", header = TRUE)

#' Basic variant filtering 
QC_filtering(data_potra_depth, data_potra_qual)
#' Filter for high-quality SNPs/Indels
snps_line26_potra <- filterVCF(vcf_potra,"results/gDNA_line26",round(lower,0),round(upper,0),c("gDNA_line3", "gDNA_line26"))

#' Extract chrom and position to a tibble
variants_line26_potra <- tibble(chrom=snps_line26_potra$`#CHROM`, pos=snps_line26_potra$POS)

#' Expand +/- 50bp around the SNPs/Indels
expanded_sequences_line26_potra <- variants_line26_potra %>%
  rowwise() %>%
  mutate(expanded_sequences = expand_sequences(chrom, pos, ref_genome))

#' View the result
print(expanded_sequences_line26_potra)

#' Retrieve just the start and stop position as a bed file for bedtools
bed_sequences_line26_potra <- tibble(expanded_sequences_line26_potra) %>% 
  mutate(start=expanded_sequences[,1], end=expanded_sequences[,2]) %>% 
  dplyr::select(chrom, start, end)

#' Write to file
write_tsv(bed_sequences_line26_potra, "results/snps_indels_line26_Potra02.bed")


#' ## Second attempt with bcftools

bcf_file <- "data/WGS/bcftoolstest/calls.bcf"

#' ## QC
calls_bcf_depth <- "data/WGS/bcftoolstest/calls_bcf_depth.txt"
data_bcf_depth <- read.delim(calls_bcf_depth, header = FALSE, sep = " ", col.names = c("CHROM", "POS", "SUM_DEPTH", "SUMSQ_DEPTH"))

system(paste("srun -A u2023008 apptainer exec -B /mnt:/mnt /mnt/picea/storage/singularity/vcftools_v0.1.16.sif vcftools --bcf", bcf_file,"--site-quality 2>/dev/null"))
data_bcf_qual <- read.table("out.lqual", header = TRUE)

QC_filtering(data_bcf_depth,data_bcf_qual)

