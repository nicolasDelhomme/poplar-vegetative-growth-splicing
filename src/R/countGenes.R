#Additional request from Daniella
#1. # of expressed genes
#2. gene ID for common/specific DE 

library(data.table)
library(DESeq2)
library(here)
library(tidyverse)

load(here("analysis/salmon/dds_mergeTechRep.rda"))

counts <- as.data.frame(counts(dds))
nrow(counts) #37075
sum(rowSums(counts) == 0) #7516
sum(rowSums(counts) > 0) #29559

sum(rowSums(counts[,grepl("T89",colnames(counts))]) > 0) #27314
sum(rowSums(counts[,grepl("line26",colnames(counts))]) > 0) #28268
sum(rowSums(counts[,grepl("line3",colnames(counts))]) > 0) #28349

sum(rowSums(counts[,grepl("T89",colnames(counts)) | grepl("line26",colnames(counts))]) > 0) #28946
sum(rowSums(counts[,grepl("T89",colnames(counts)) | grepl("line3",colnames(counts))]) > 0) #28955
sum(rowSums(counts[,grepl("T89",colnames(counts))]) > 0 & rowSums(counts[,grepl("line26",colnames(counts))]) > 0) #26636
sum(rowSums(counts[,grepl("T89",colnames(counts))]) > 0 & rowSums(counts[,grepl("line3",colnames(counts))]) > 0) #26708
length(rownames(counts)[rowSums(counts[,grepl("T89",colnames(counts))]) == 0 & rowSums(counts[,grepl("line3",colnames(counts))]) > 0]) #1641
sum(rownames(counts)[rowSums(counts[,grepl("T89",colnames(counts))]) == 0 & rowSums(counts[,grepl("line3",colnames(counts))]) > 0] %in% l3$...1) #205
length(rownames(counts)[rowSums(counts[,grepl("T89",colnames(counts))]) == 0 & rowSums(counts[,grepl("line26",colnames(counts))]) > 0]) #1632
sum(rownames(counts)[rowSums(counts[,grepl("T89",colnames(counts))]) == 0 & rowSums(counts[,grepl("line26",colnames(counts))]) > 0] %in% l26$...1) #204

l26 <- read_csv(here("analysis/DE/Line26_genes.csv"))
l3 <- read_csv(here("analysis/DE/Line03_genes.csv"))
nrow(l26) #8765
nrow(l3) #8064

length(intersect(l26$...1,l3$...1)) #5292
write.table(intersect(l26$...1,l3$...1),here("analysis/DE/commonL26L3_all.txt"),quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(intersect(l26[l26$log2FoldChange > 0,]$...1,l3[l3$log2FoldChange > 0,]$...1),here("analysis/DE/commonL26L3_up.txt"),quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(intersect(l26[l26$log2FoldChange < 0,]$...1,l3[l3$log2FoldChange < 0,]$...1),here("analysis/DE/commonL26L3_down.txt"),quote = FALSE, row.names = FALSE, col.names = FALSE)

length(setdiff(l26$...1,l3$...1)) #3473
write.table(setdiff(l26$...1,l3$...1),here("analysis/DE/specificL26_all.txt"),quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(setdiff(l26[l26$log2FoldChange > 0,]$...1,l3[l3$log2FoldChange > 0,]$...1),here("analysis/DE/specificL26_up.txt"),quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(setdiff(l26[l26$log2FoldChange < 0,]$...1,l3[l3$log2FoldChange < 0,]$...1),here("analysis/DE/specificL26_down.txt"),quote = FALSE, row.names = FALSE, col.names = FALSE)

length(setdiff(l3$...1,l26$...1)) #2772
write.table(setdiff(l3$...1,l26$...1),here("analysis/DE/specificL3_all.txt"),quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(setdiff(l3[l3$log2FoldChange > 0,]$...1,l26[l26$log2FoldChange > 0,]$...1),here("analysis/DE/specificL3_up.txt"),quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(setdiff(l3[l3$log2FoldChange < 0,]$...1,l26[l26$log2FoldChange < 0,]$...1),here("analysis/DE/specificL3_down.txt"),quote = FALSE, row.names = FALSE, col.names = FALSE)
