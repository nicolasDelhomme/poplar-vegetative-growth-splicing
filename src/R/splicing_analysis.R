suppressPackageStartupMessages({
  library(ASpli)
  library(GenomicFeatures)
  library(readxl)
  library(openxlsx)
  library(tidyverse)
  library(here)
  library(gplots)
  library(RColorBrewer)
})

suppressMessages({
  source(here("UPSCb-common/src/R/topGoUtilities.R"))
  source(here("UPSCb-common/src/R/featureSelection.R"))
})


#' * Functions
#' 

#' Extract results and make GO enrichment plots
extractEnrichmentResults <- function(enrichment,
                                     go.namespace=c("BP","CC","MF"),
                                     genes=NULL,export=FALSE,plot=TRUE,
                                     default_dir=here("analysis/splicing_analysis"),
                                     default_prefix="DAS-",
                                     sample_name=NULL){
  # sanity
  if(is.null(unlist(enrichment)) | length(unlist(enrichment)) == 0){
    message("No GO enrichment for",names(enrichment))
  } else {
    # write out
    if(export){
      write_tsv(bind_rows(enrichment),
                file=here(default_dir,
                          paste0(default_prefix,sample_name,"-genes_GO-enrichment.tsv")))
      if(!is.null(genes)){
        write_tsv(
          enrichedTermToGenes(genes=genes,terms=enrichment$id,url=url,mc.cores=16L),
          file=here(default_dir,
                    paste0(default_prefix,sample_name,"-enriched-term-to-genes.tsv"))
        )
      }
    }
    if(plot){
      gocatname <- c(BP="Biological Process",
                     CC="Cellular Component",
                     MF="Molecular Function")
      
      for (n in go.namespace) {
        dat <- enrichment[[n]]
        # Print dat to check its content
        print(dat)
        
        # Rest of your code for plotting
        dat$GeneRatio <- dat$Significant/dat$Annotated
        dat$adjustedPvalue <- as.numeric(dat$FDR)
        dat$Count <- as.numeric(dat$Significant)
        dat <- dat[order(dat$GeneRatio),]
        dat$Term <- factor(dat$Term, levels = unique(dat$Term))
        
        print(
          ggplot(dat, aes(x = Term, y = GeneRatio, color = adjustedPvalue, size = Count)) + 
            geom_point() +
            scale_color_gradient(low = "red", high = "blue") +
            theme_bw() + 
            ylab("GeneRatio") + 
            xlab("") + 
            ggtitle(paste0("GO enrichment: ", gocatname[n], " for sample: ", sample_name)) +
            coord_flip()
        )
      }
    }
  }
}

#'  Define a function to process and export data for each GO term
process_GO_term <- function(term_type) {
  # Create data frame for the current GO term
  export_enr <- data.frame(
    GO.ID = term_type$GO.ID, 
    Term = term_type$Term, 
    Annotated = term_type$Annotated,
    Significant = term_type$Significant,
    Expected = term_type$Expected,
    FDR = term_type$FDR,
    Significant_genes = sapply(term_type$siggenes, function(x) paste(unlist(x), collapse = "|"))
  )
  return(export_enr)
}

#' # Alternative splicing analysis
#' 
#' 

#' GFF preprocessing
gtfFileName <- "reference/gff/Potra02_genes.gff"
genomeTxDb <- makeTxDbFromGFF( gtfFileName )

saveDb(genomeTxDb,file="Potra_genes.sqlite")

#' Feature extraction
features <- binGenome( genomeTxDb )

#' Prepare BAMS and target file
bams <- list.files("/mnt/ada/projects/aspseq/mschmid/poplar-CRISPR-WGS/RNA-seq/preprocessed/STAR", pattern="*.bam$", full.names = TRUE)

targets <- data.frame(row.names = paste0(rep(c("L03", "L26", "T89"), each = 4), rep("_", times = 12), rep(1:4, times = 3)),
                      bam = bams,
                      genotype = c( 'line3','line3','line3','line3','line26','line26','line26','line26','control','control','control','control'),
                      stringsAsFactors = FALSE)

mBAMs <- data.frame( bam = list.files("/mnt/ada/projects/aspseq/mschmid/poplar-CRISPR-WGS/RNA-seq/preprocessed/STAR/merged", pattern="*.bam$", full.names = TRUE),
                     condition = c("line3","line26", "control"))

#' Read counting against annotated features
gbcounts <- gbCounts(features=features, targets=targets,
                     minReadLength = 100, maxISize = 50000, libType = "PE", strandMode = 2)
gbcounts

#' Junction-based de novo counting and splicing signal estimation
asd <- jCounts(counts=gbcounts, features=features, minReadLength=100)

asd

#' Differential splicing expression
gb_line3 <- gbDUreport(gbcounts, contrast = c(1,0,-1))
gb_line26 <- gbDUreport(gbcounts, contrast = c(0,1,-1))

#' Differential junction usage analysis
jdur_line3 <- jDUreport(asd, contrast=c(1,0,-1))
jdur_line26 <- jDUreport(asd, contrast=c(0,1,-1))

#' Bin and junction signal integration
sr_line3 <- splicingReport(gb_line3, jdur_line3, counts=gbcounts)
sr_line26 <- splicingReport(gb_line26, jdur_line26, counts=gbcounts)
sr <- splicingReport(gb_t89, jdur_t89, counts=gbcounts)

#' Summary of integration of splicing signals along genomic-regions
is_line3 <- integrateSignals(sr_line3,asd)
is_line26 <- integrateSignals(sr_line26,asd)

#' Export results
exportIntegratedSignals(is_line3,sr=sr_line3,
                        output.dir = "Line3_AS_stranded",
                        counts=gbcounts,features=features,asd=asd,
                        mergedBams = mBAMs)

exportIntegratedSignals(is_line26,sr=sr_line26,
                        output.dir = "Line26_AS_stranded",
                        counts=gbcounts,features=features,asd=asd,
                        mergedBams = mBAMs)


#' # Summarize splicing events per mutant
#' 
#' ## Line 3
#' 
#' #' Read in Excel
as_line3 <- read_excel("doc/line3_control_stranded.xlsx")
#' Get the frequency per splicing event
events_line3 <- as.data.frame(table(as_line3$Event))
#' Remove unassigned
events_line3 <- events_line3[-1,]
#' Get total number of splicing events
total_events_line3 <- sum(events_line3$Freq)
#' Calculate percentages
events_line3$Percentage <- (events_line3$Freq / total_events_line3) * 100
#' Combine annotations
events_line3$Annotation <- c("Alt 3'", "Alt 3'", "Alt 5'", "Alt 5'", "ASCE", "ES", "ES", "IoR", "IR", "IR", "Novel", "Novel", "Novel", "Novel", "Undefined")
#' Add annotation for line for plotting
events_line3$Line <- rep("Line3",15)
#' Filter some types of events
events_line3 <- events_line3[which(!events_line3$Annotation %in% c("ASCE", "Undefined", "IoR", "Novel")),]

#' ## Line 26
#' 
#' Read in Excel
as_line26 <- read_excel("doc/line26_control_stranded.xlsx")
#' Get the frequency per splicing event
events_line26 <- as.data.frame(table(as_line26$Event))
#' Remove unassigned
events_line26 <- events_line26[-1,]
#' Get total number of splicing events
total_events_line26 <- sum(events_line26$Freq)
#' Calculate percentages
events_line26$Percentage <- (events_line26$Freq / total_events_line26) * 100
#' Combine annotations
events_line26$Annotation <- c("Alt 5'/3'","Alt 3'", "Alt 3'", "Alt 5'", "Alt 5'", "ASCE", "ES", "ES", "IoR", "IR", "IR", "Novel", "Novel", "Novel", "Novel", "Undefined")
#' Add annotation for line for plotting
events_line26$Line <- rep("Line26",16)
#' Filter some types of events
events_line26 <- events_line26[which(!events_line26$Annotation %in% c("ASCE", "Undefined", "IoR", "Novel")),]

#' ## Get the total annotated events¨
#' 

#' Local AS events annotated in the Potra genome (from SUPPA2)

es_anno <- 2874
ir_anno <- 7215
three_anno <- 4128
five_anno <- 5080

total_as_annotated <- es_anno + ir_anno + three_anno + five_anno

events_as_anno <- data.frame(
  Var1 = c("ES","IR", "Alt 3'", "Alt 5'"),
  Freq = c(es_anno, ir_anno, three_anno, five_anno),
  Annotation = c("ES", "IR", "Alt 3'", "Alt 5'"),
  Percentage = c(es_anno, ir_anno, three_anno, five_anno) / total_as_annotated, 
  Line = rep("Total Annotated", 4)
)


#' Combine events for Line 3 and Total Annotated
events_line3_anno <- rbind(events_line3,events_as_anno)

#' Plot 
ggplot(data = events_line3_anno, aes(x = Line, y = Percentage, fill = Annotation,
                              label = sprintf("%.02f", Percentage))) + geom_bar(position = "fill", stat = "identity") +
  scale_fill_viridis_d() +  # Colorblind-friendly and print-friendly palette
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  labs(x = "", y = "% AS events")

#' Combine events for Line 3 and Total Annotated
events_line26_anno <- rbind(events_line26,events_as_anno)

#' Plot 
ggplot(data = events_line26_anno, aes(x = Line, y = Percentage, fill = Annotation,
                                     label = sprintf("%.02f", Percentage))) + geom_bar(position = "fill", stat = "identity") +
  scale_fill_viridis_d() +  # Colorblind-friendly and print-friendly palette
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  labs(x = "", y = "% AS events")


#' ## Extract the unique DAS genes
#' 
DAS_genes_line3 <- unique(as_line3$Locus)

DAS_genes_line26 <- unique(as_line26$Locus)

write_csv(as.data.frame(DAS_genes_line3), "DAS_genes_line3.csv")

print(paste0("Number of AS regulated genes: ", length(DAS_genes_line3)))

write_csv(as.data.frame(DAS_genes_line26), "DAS_genes_line26.csv")

print(paste0("Number of AS regulated genes: ", length(DAS_genes_line26)))


#' ## Compare with DE genes
#' 
#' Line 3

DEG_line3 <- read_csv("analysis/DE/Line03_genes.csv") 

colnames(DEG_line3)[1] <- "Locus"

common_genes_line3 <- intersect(as_line3$Locus, DEG_line3$Locus)

print(paste0("Number of Transcription and AS regulated genes: ", length(common_genes_line3)))

write_csv(as.data.frame(common_genes_line3), "DE_AS_genes_line3.csv")
#' 
#' Line 26

DEG_line26 <- read_csv("analysis/DE/Line26_genes.csv") 

colnames(DEG_line26)[1] <- "Locus"

common_genes_line26 <- intersect(as_line26$Locus, DEG_line26$Locus)

print(paste0("Number of Transcription and AS regulated genes: ", length(common_genes_line26)))

write_csv(as.data.frame(common_genes_line26), "DE_AS_genes_line26.csv")


#' # GO enrichment using TopGO
#' 

goannot <- prepAnnot(mapping = "/mnt/picea/storage/reference/Populus-tremula/v2.2/gopher/gene_to_go.tsv")
p_val_cutoff= 0.05
padj_method = "none"

#' ## Line3
#' 
#background_line3 <- setNames(gb_line3@genes$pvalue, rownames(gb_line3@genes))
background_line3 <- rownames(gb_line3@genes)
enr_line3 <- topGO(DAS_genes_line3,background_line3,goannot,alpha=p_val_cutoff,p.adjust=padj_method,getgenes=TRUE)

extractEnrichmentResults(enr_line3, sample_name="Line3")

#' List of GO terms (BP, CC, MF)
go_terms_line3 <- list(BP = enr_line3$BP, CC = enr_line3$CC, MF = enr_line3$MF)

exported_data_line3 <- lapply(go_terms_line3, process_GO_term)

# Create a new Excel workbook
wb_line3_1 <- createWorkbook()

# Apply the process_GO_term function to each GO term
lapply(names(exported_data_line3), function(go_term) {
  # Create data frame for the current GO term
  export_line3 <- exported_data_line3[[go_term]]
  
  # Add the data frame to a new sheet in the Excel workbook
  addWorksheet(wb_line3_1, sheetName = go_term)
  writeData(wb_line3_1, sheet = go_term, x = export_line3)
})

# Save the Excel workbook
saveWorkbook(wb_line3_1, "GO_analysis_DAS_genes_Line3.xlsx", overwrite = TRUE)


wb_line3_2 <- createWorkbook()

# Iterate over each list element
for (term in names(exported_data_line3)) {
  
  # Extract data for the current term
  term_data_line3 <- exported_data_line3[[term]]
  
  # Perform separation of genes and arrange the data
  term_data_separated_line3 <- tidyr::separate_rows(term_data_line3, Significant_genes, sep = "\\|")
  term_data_sorted_line3 <- term_data_separated_line3 %>%
    group_by(Significant_genes) %>%
    arrange(Significant_genes) %>%
    summarize(GO.ID = paste(GO.ID, collapse = "|"))
  
  # Add the sorted data to a new sheet in the Excel workbook
  addWorksheet(wb_line3_2, sheetName = term)
  writeData(wb_line3_2, sheet = term, x = term_data_sorted_line3)
}

# Save the Excel workbook
saveWorkbook(wb_line3_2, "GO_terms_per_gene_Line3.xlsx", overwrite = TRUE)

#' ## Line26
#' 
background_line26 <- rownames(gb_line26@genes)

enr_line26 <- topGO(DAS_genes_line26,background_line26,goannot,alpha=p_val_cutoff,p.adjust=padj_method,getgenes=TRUE)

extractEnrichmentResults(enr_line26, sample_name="Line26")

go_terms_line26 <- list(BP = enr_line26$BP, CC = enr_line26$CC, MF = enr_line26$MF)

# Apply the process_GO_term function to each GO term
exported_data_line26 <- lapply(go_terms_line26, process_GO_term)

# Create a new Excel workbook
wb_line26_1 <- createWorkbook()

# Apply the process_GO_term function to each GO term
lapply(names(exported_data_line26), function(go_term) {
  # Create data frame for the current GO term
  export_enr <- exported_data_line26[[go_term]]
  
  # Add the data frame to a new sheet in the Excel workbook
  addWorksheet(wb_line26_1, sheetName = go_term)
  writeData(wb_line26_1, sheet = go_term, x = export_enr)
})

# Save the Excel workbook
saveWorkbook(wb_line26_1, "GO_analysis_DAS_genes_Line26.xlsx", overwrite = TRUE)

wb_line26_2 <- createWorkbook()

# Iterate over each list element
for (term in names(exported_data_line26)) {
  
  # Extract data for the current term
  term_data_line26 <- exported_data_line26[[term]]
  
  # Perform separation of genes and arrange the data
  term_data_separated_line26 <- tidyr::separate_rows(term_data_line26, Significant_genes, sep = "\\|")
  term_data_sorted_line26 <- term_data_separated_line26 %>%
    group_by(Significant_genes) %>%
    arrange(Significant_genes) %>%
    summarize(GO.ID = paste(GO.ID, collapse = "|"))
  
  # Add the sorted data to a new sheet in the Excel workbook
  addWorksheet(wb_line26_2, sheetName = term)
  writeData(wb_line26_2, sheet = term, x = term_data_sorted_line26)
}

# Save the Excel workbook
saveWorkbook(wb_line26_2, "GO_terms_per_gene_Line26.xlsx", overwrite = TRUE)
