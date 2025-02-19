#' ---
#' title: "Poplar CRISPR lines Differential Expression"
#' author: "Fai Theerarat Kochakarn"
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
    library(DESeq2)
    library(gplots)
    library(here)
    library(hyperSpec)
    library(RColorBrewer)
    library(tidyverse)
    library(VennDiagram)
    library(emoji)
    library(ggVennDiagram)
})

#' * Helper files
suppressMessages({
    source(here("UPSCb-common/Rtoolbox/src/plotEnrichedTreemap.R"))
    source(here("UPSCb-common/src/R/featureSelection.R"))
    source(here("UPSCb-common/src/R/volcanoPlot.R"))
    source(here("UPSCb-common/src/R/gopher.R"))
})

#' * Graphics
pal=brewer.pal(8,"Dark2")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

#' * Functions
#' 1. plot specific gene expression
#' ```{r edit1, echo=FALSE,eval=FALSE}
#' CHANGEME - here you need to change the variables in the 
#' plot to display the expression values accross your samples
#' The example below has 2 variables MGenotype and MDay. These 
#' need replacing by the variable(s) of interest in your project
#' ```
"line_plot" <- function(dds=dds,vst=vst,gene_id=gene_id){
    message(paste("Plotting",gene_id))
    sel <- grepl(gene_id,rownames(vst))
    stopifnot(sum(sel)==1)

    p <- ggplot(bind_cols(as.data.frame(colData(dds)),
                          data.frame(value=vst[sel,])),
                aes(x=Genotype,y=value,fill=Genotype)) +
        geom_boxplot() +
        geom_jitter(color="black") +
        scale_y_continuous(name="VST expression") + 
        ggtitle(label=paste("Expression for: ",gene_id))
    
    suppressMessages(suppressWarnings(plot(p)))
    return(NULL)
}

#' 2. extract the DE results. Default cutoffs are
#' from Schurch _et al._, RNA, 2016
"extract_results" <- function(dds,vst,contrast,
                              padj=0.01,lfc=0.5,
                              plot=TRUE,verbose=TRUE,
                              export=TRUE,default_dir=here("data/analysis/DE"),
                              default_prefix="DE-",
                              labels=colnames(dds),
                              sample_sel=1:ncol(dds),
                              expression_cutoff=0,
                              debug=FALSE,filter=c("median",NULL),...){
    
    # get the filter
    if(!is.null(match.arg(filter))){
        filter <- rowMedians(counts(dds,normalized=TRUE))
        message("Using the median normalized counts as default, set filter=NULL to revert to using the mean")
    }
    
    # validation
    if(length(contrast)==1){
        res <- results(dds,name=contrast,filter = filter)
    } else {
        res <- results(dds,contrast=contrast,filter = filter)
    }
    
    stopifnot(length(sample_sel)==ncol(vst))
    
    if(plot){
        par(mar=c(5,5,5,5))
        volcanoPlot(res)
        par(mar=mar)
    }
    
    # a look at independent filtering
    if(plot){
        plot(metadata(res)$filterNumRej,
             type="b", ylab="number of rejections",
             xlab="quantiles of filter")
        lines(metadata(res)$lo.fit, col="red")
        abline(v=metadata(res)$filterTheta)
    }
    
    if(verbose){
        message(sprintf("The independent filtering cutoff is %s, removing %s of the data",
                        round(metadata(res)$filterThreshold,digits=5),
                        names(metadata(res)$filterThreshold)))
        
        max.theta <- metadata(res)$filterNumRej[which.max(metadata(res)$filterNumRej$numRej),"theta"]
        message(sprintf("The independent filtering maximises for %s %% of the data, corresponding to a base mean expression of %s (library-size normalised read)",
                        round(max.theta*100,digits=5),
                        round(quantile(counts(dds,normalized=TRUE),probs=max.theta),digits=5)))
    }
    
    if(plot){
        qtl.exp=quantile(counts(dds,normalized=TRUE),probs=metadata(res)$filterNumRej$theta)
        dat <- data.frame(thetas=metadata(res)$filterNumRej$theta,
                          qtl.exp=qtl.exp,
                          number.degs=sapply(lapply(qtl.exp,function(qe){
                              res$padj <= padj & abs(res$log2FoldChange) >= lfc & 
                                  ! is.na(res$padj) & res$baseMean >= qe
                          }),sum))
        if(debug){
            plot(ggplot(dat,aes(x=thetas,y=qtl.exp)) + 
                     geom_line() + geom_point() +
                     scale_x_continuous("quantiles of expression") + 
                     scale_y_continuous("base mean expression") +
                     geom_hline(yintercept=expression_cutoff,
                                linetype="dotted",col="red"))
        
            p <- ggplot(dat,aes(x=thetas,y=qtl.exp)) + 
                geom_line() + geom_point() +
                scale_x_continuous("quantiles of expression") + 
                scale_y_log10("base mean expression") + 
                geom_hline(yintercept=expression_cutoff,
                           linetype="dotted",col="red")
            suppressMessages(suppressWarnings(plot(p)))
            
            plot(ggplot(dat,aes(x=thetas,y=number.degs)) + 
                     geom_line() + geom_point() +
                     geom_hline(yintercept=dat$number.degs[1],linetype="dashed") +
                     scale_x_continuous("quantiles of expression") + 
                     scale_y_continuous("Number of DE genes"))
            
            plot(ggplot(dat,aes(x=thetas,y=number.degs[1] - number.degs),aes()) + 
                     geom_line() + geom_point() +
                     scale_x_continuous("quantiles of expression") + 
                     scale_y_continuous("Cumulative number of DE genes"))
            
            plot(ggplot(data.frame(x=dat$thetas[-1],
                                   y=diff(dat$number.degs[1] - dat$number.degs)),aes(x,y)) + 
                     geom_line() + geom_point() +
                     scale_x_continuous("quantiles of expression") + 
                     scale_y_continuous("Number of DE genes per interval"))
            
            plot(ggplot(data.frame(x=dat$qtl.exp[-1],
                                   y=diff(dat$number.degs[1] - dat$number.degs)),aes(x,y)) + 
                     geom_line() + geom_point() +
                     scale_x_continuous("base mean of expression") + 
                     scale_y_continuous("Number of DE genes per interval"))
            
            p <- ggplot(data.frame(x=dat$qtl.exp[-1],
                                   y=diff(dat$number.degs[1] - dat$number.degs)),aes(x,y)) + 
                geom_line() + geom_point() +
                scale_x_log10("base mean of expression") + 
                scale_y_continuous("Number of DE genes per interval") + 
                geom_vline(xintercept=expression_cutoff,
                           linetype="dotted",col="red")
            suppressMessages(suppressWarnings(plot(p)))
        }
    }
    
    sel <- res$padj <= padj & abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & 
        res$baseMean >= expression_cutoff
    
    if(verbose){
      message(sprintf(paste(
        ifelse(sum(sel)==1,
               "There is %s gene that is DE",
               "There are %s genes that are DE"),
        "with the following parameters: FDR <= %s, |log2FC| >= %s, base mean expression > %s"),
        sum(sel),padj,
        lfc,expression_cutoff))
    }
    
    # proceed only if there are DE genes
    if(sum(sel) > 0){
        val <- rowSums(vst[sel,sample_sel,drop=FALSE])==0
        if (sum(val) >0){
          warning(sprintf(paste(
            ifelse(sum(val)==1,
                   "There is %s DE gene that has",
                   "There are %s DE genes that have"),
            "no vst expression in the selected samples"),sum(val)))
          sel[sel][val] <- FALSE
        } 

        if(export){
            if(!dir.exists(default_dir)){
                dir.create(default_dir,showWarnings=FALSE,recursive=TRUE,mode="0771")
            }
            write.csv(res,file=file.path(default_dir,paste0(default_prefix,"results.csv")))
            write.csv(res[sel,],file.path(default_dir,paste0(default_prefix,"genes.csv")))
        }
        if(plot & sum(sel)>1){
            heatmap.2(t(scale(t(vst[sel,sample_sel]))),
                      distfun = pearson.dist,
                      hclustfun = function(X){hclust(X,method="ward.D2")},
                      trace="none",col=hpal,labRow = FALSE,
                      labCol=labels[sample_sel],...
            )
        }
    }
    return(list(all=rownames(res[sel,]),
                up=rownames(res[sel & res$log2FoldChange > 0,]),
                dn=rownames(res[sel & res$log2FoldChange < 0,])))
}

#' 3. extract and plot the enrichment results
extractEnrichmentResults <- function(enrichment,task="go",
                                     diff.exp=c("all","up","dn"),
                                     go.namespace=c("BP","CC","MF"),
                                     genes=NULL,export=TRUE,plot=TRUE,
                                     default_dir=here("data/analysis/DE"),
                                     default_prefix="DE",
                                     url="athaliana"){
    # process args
    diff.exp <- match.arg(diff.exp)
    de <- ifelse(diff.exp=="all","none",
                 ifelse(diff.exp=="dn","down",diff.exp))

    # sanity
    if( is.null(enrichment[[task]]) | length(enrichment[[task]]) == 0){
        message(paste("No enrichment for",task))
    } else {

        # write out
        if(export){
            write_tsv(enrichment[[task]],
                      file=here(default_dir,
                                paste0(default_prefix,"-genes_GO-enrichment.tsv")))
            if(!is.null(genes)){
                write_tsv(
                    enrichedTermToGenes(genes=genes,terms=enrichment[[task]]$id,url=url,mc.cores=16L),
                    file=here(default_dir,
                              paste0(default_prefix,"-enriched-term-to-genes.tsv"))
                )
            }
        }
        
        if(plot){
            sapply(go.namespace,function(ns){
                titles <- c(BP="Biological Process",
                            CC="Cellular Component",
                            MF="Molecular Function")
                suppressWarnings(tryCatch({plotEnrichedTreemap(enrichment,enrichment=task,
                                                               namespace=ns,
                                                               de=de,title=paste(default_prefix,titles[ns]))},
                                          error = function(e) {
                                              message(paste("Treemap plot failed for",ns, 
                                                            "because of:",e))
                                              return(NULL)
                                          }))
            })
        }
    }
}

#' * Data
#' ```{r load, echo=FALSE,eval=FALSE}
#' CHANGEME - here you are meant to load an RData object
#' that contains a DESeqDataSet object. If you ran the 
#' biological QA template, you need not change anything
#' ```
load(here("analysis/salmon/dds_mergeTechRep.rda"))
dds$Genotype <- relevel(dds$Genotype,ref = "T89")

#' ## Normalisation for visualisation
vsd <- varianceStabilizingTransformation(dds,blind=FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)
dir.create(here("analysis/DE"),showWarnings=FALSE)
save(vst,file=here("analysis/DE/vst-aware.rda"))
write_delim(as.data.frame(vst) %>% rownames_to_column("ID"),
            here("analysis/DE/vst-aware.tsv"))

#' ## Gene of interests
#' ```{r goi, echo=FALSE,eval=FALSE}
#' CHANGEME - Here, you can plot the expression pattern of your gene of
#' interest. You need to have a list of genes in a text file, one geneID per line
#' The ID should exist in your vst data.
#' Note that for the plot to work, you also need to edit the first function (`line_plot`)
#' at the top of this file
#' ```
#' ### Knocked-out gene
#' The targeted locus for line 26 is `Potri.006G174000/Potra001877g14982` (Potra2n6c13821)
#' And for the line 03, it is `Potri.018G096200/Potra001273g10998` (Potra2n6c13821)
#' This is the expression level of targeted gene.
kogene <- "Potra2n6c13821"
stopifnot(kogene %in% rownames(vst))
line_plot(dds=dds,vst=vst,gene_id = kogene)

#' ### Other SM and LSM
suppressMessages(goi <- read_delim(here("doc/goi.txt")))
knitr::kable(goi)
genelist <- goi$V2[!(goi$V2 == "Potra2n6c13821")]
stopifnot(all(genelist %in% rownames(vst)))
dev.null <- lapply(genelist,line_plot,dds=dds,vst=vst)

#' ## Differential Expression
dds <- DESeq(dds)

#' * Dispersion estimation
#' The dispersion estimation is adequate
plotDispEsts(dds)

#' Check the different contrasts
resultsNames(dds)

#' ## Results
#' ```{r res, echo=FALSE,eval=FALSE}
#' CHANGEME - here you need to define the contrast you want to 
#' study - see the example in the next block. 
#' 
#' The `contrast` can be given
#' by name, as a list (numerator/denominator) or as a vector of weight (e.g. c(0,1));
#' read the DESeq2 vignette for more info
#' 
#' The `label` argument is typically one (or a combination) of the metadata stored in colData
#' 
#' The function allows for looking at the independent filtering results using `debug=TRUE`
#' 
#' If you are not satisfied with the default from DESeq2, you can set your own cutoff using `expression_cutoff`
#' 
#' You can change the default output file prefix using `default_prefix`
#' 
#' You can select the set of samples to be added to the `heatmap`, using the `sample_sel` argument. It takes a logical vector.
#' 
#' ```

#' #### Line 26
line26 <- extract_results(dds=dds,vst=vst,contrast="Genotype_line26_vs_T89",
														default_dir = here("analysis/DE"),
														default_prefix = "Line26_",
														labels=dds$BioID,
														sample_sel = dds$Genotype %in% c("line26","T89"), 
														cexCol=.9, plot = TRUE, verbose = TRUE)

#' #### Line 03
line03 <- extract_results(dds=dds,vst=vst,contrast="Genotype_line3_vs_T89",
														 default_dir = here("analysis/DE"),
														 default_prefix = "Line03_",
														 labels=dds$BioID,
														 sample_sel = dds$Genotype %in% c("line3","T89"), 
														 cexCol=.9, plot = TRUE, verbose = TRUE)

#' ### Venn Diagram
#' ```{r venn, echo=FALSE,eval=FALSE}
#' CHANGEME - Here, you typically would have run several contrasts and you want
#' to assess their overlap plotting VennDiagrams.
#' 
#' In the examples below, we assume that these resutls have been saved in a list
#' called `res.list`
#' ```

res.list <- list("Line 26"=line26,
								 "Line 03"=line03)

#' #### All DE genes
ggVennDiagram(lapply(res.list,"[[","all"), label_alpha = 0, label = "count") + 
	scale_fill_gradient(low="white",high = "dodgerblue") + 
	scale_x_continuous(expand = expansion(mult = .1))

#' #### Up-regulated DE genes (up in treated)
ggVennDiagram(lapply(res.list,"[[","up"), label_alpha = 0, label = "count") + 
	scale_fill_gradient(low="white",high = "dodgerblue") + 
	scale_x_continuous(expand = expansion(mult = .1))

#' #### Down-regulated DE genes (down in treated)
ggVennDiagram(lapply(res.list,"[[","dn"), label_alpha = 0, label = "count") + 
	scale_fill_gradient(low="white",high = "dodgerblue") + 
	scale_x_continuous(expand = expansion(mult = .1))


#' ### Gene Ontology enrichment
#' ```{r go, echo=FALSE,eval=FALSE}
#' Once you have obtained a list of candidate genes, you most probably want
#' to annotate them.
#' 
#' In the following example, we first identify the background; _i.e._ the
#' population of expressed genes. We select the genes expressed in a least
#' 2 replicate of one condition at a cutoff of `exp`.
#' 
#' Next we run the enrichment, in the example against `athaliana` using 
#' the gofer3 REST API (interfaced through the gopher.R script loaded at the
#' beginning of this fil).
#' 
#' Finally we export the go enrichment as a complete table.
#' We used to export another table consisting
#' of only the `id` and `padj` columns for using as input for _e.g._
#' REVIGO; but since flash is EOL and REVIGO not updated, we instead rely on 
#' the RtoolBox treemap.
#' 
#' In addition we now also export the list of genes that most likely resulted in
#' the corresponding go enrichment.
#' 
#' Make sure to change the `url` to match your species
#' 
#' ```
background <- rownames(vst)[featureSelect(vst,dds$Genotype,exp=0.00001)]

#' #### DEGs from Line26

res.list <- list("Line 26"=line26)

enr.list <- lapply(res.list,function(r){
	lapply(r,gopher,background=background,task="go",url="potra2")
})

dev.null <- lapply(names(enr.list),function(n){
	lapply(names(enr.list[[n]]),function(de){
		extractEnrichmentResults(enr.list[[n]][[de]],
														 diff.exp=de,
														 genes=res.list[[n]][[de]],
														 default_prefix=paste(n,de,sep="-"),
														 url="potra2")
	})
})

#' #### DEGs from Line03

res.list <- list("Line 03"=line03)

enr.list <- lapply(res.list,function(r){
	lapply(r,gopher,background=background,task="go",url="potra2")
})

dev.null <- lapply(names(enr.list),function(n){
	lapply(names(enr.list[[n]]),function(de){
		extractEnrichmentResults(enr.list[[n]][[de]],
														 diff.exp=de,
														 genes=res.list[[n]][[de]],
														 default_prefix=paste(n,de,sep="-"),
														 url="potra2")
	})
})

#' #### Common DEGs between two lines
res.list <- list("CommonLine26andLine03"=list("all"=intersect(line26$all,line03$all),
																							"up"=intersect(line26$up,line03$up),
																							"dn"=intersect(line26$dn,line03$dn)))

enr.list <- lapply(res.list,function(r){
	lapply(r,gopher,background=background,task="go",url="potra2")
})

dev.null <- lapply(names(enr.list),function(n){
	lapply(names(enr.list[[n]]),function(de){
		extractEnrichmentResults(enr.list[[n]][[de]],
														 diff.exp=de,
														 genes=res.list[[n]][[de]],
														 default_prefix=paste(n,de,sep="-"),
														 url="potra2")
	})
})

#' # Session Info 
#'  ```{r session info, echo=FALSE}
#'  sessionInfo()
#'  ```


