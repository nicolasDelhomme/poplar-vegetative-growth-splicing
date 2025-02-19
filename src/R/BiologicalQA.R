#' ---
#' title: "CRISPR lines Biological QA"
#' subtitle: "Porcupine-like lines"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    fig_width: 9
#'    fig_height: 6
#'    toc: true
#'    number_sections: true
#'    toc_depth: 3
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
#' This section and the next are relevant for reproducibility purposes. For results, please skip to section 3 (Quality Control)
#' 
#' * Libraries
suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(emoji)
  library(gplots)
  library(gtools)
  library(here)
  library(hyperSpec)
  library(limma)
  library(magrittr)
  library(parallel)
  library(patchwork)
  library(PCAtools)
  library(pheatmap)
  library(plotly)
  library(pvclust)
  library(RColorBrewer)
  library(tidyverse)
  library(tximport)
  library(vsn)
})

#' * Helper functions
source(here("UPSCb-common/src/R/featureSelection.R"))

#' * Graphics
hpal <- colorRampPalette(c("blue","white","red"))(100)
pal <- brewer.pal(9,"Blues")

#' # Data
#' * Sample information
samples <- read_tsv(here("doc/samples.tsv"),
                    col_types=cols(.default=col_factor()))

#' * tx2gene translation table
tx2gene <- suppressMessages(read_delim(here("reference/annotation/tx2gene.tsv.gz"),delim="\t",
                                 col_names=c("TXID","GENE")))

#' * Raw data
filelist <- list.files(here("analysis/salmon"),
                          recursive = TRUE, 
                          pattern = "quant.sf",
                          full.names = TRUE)

#' * Sanity check to ensure that the data is sorted according to the sample info
stopifnot(all(match(sub("_1_sort.*","",basename(dirname(filelist))),
                    samples$SampleID) == 1:nrow(samples)))

#' * add filelist to samples as a new column
samples %<>% mutate(Filenames = filelist)

#' Read the expression at the gene level
txi <- suppressMessages(tximport(files = samples$Filenames,
                                 type = "salmon",
                                 tx2gene=tx2gene))
counts <- txi$counts
colnames(counts) <- samples$SampleID
#' 
#' <hr />
#' &nbsp;
#' 
#' # Quality Control
#' * "Not expressed" genes
sel <- rowSums(counts) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(counts),digits=1),
        sum(sel),
        nrow(counts))

#' ## Sequencing depth
#' * Let us take a look at the sequencing depth, colouring by CHANGEME
dat <- tibble(x=colnames(counts),y=colSums(counts)) %>% 
  bind_cols(samples)

ggplot(dat,aes(x,y,fill=Genotype)) + 
  geom_col() + 
  scale_y_continuous(name="reads") +
  facet_grid(~ factor(Genotype), scales = "free") +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=90,size=4),
        axis.title.x=element_blank())

#' `r emoji("point_right")` **We observe a lot of difference in the raw sequencing depth, but we know some are tech replicates.**
#'
#' ## per-gene mean expression
#' 
#' _i.e._ the mean raw count of every gene across samples is calculated
#' and displayed on a log10 scale.
#' 

ggplot(data.frame(value=log10(rowMeans(counts))),aes(x=value)) + 
  geom_density(na.rm = TRUE) +
  ggtitle("gene mean raw counts distribution") +
  scale_x_continuous(name="mean raw counts (log10)") + 
  theme_bw()

#' `r emoji("point_right")` **The cumulative gene coverage is as expected, despite the lower coverage of some samples**
#' 
#' ## Per-sample expression

dat <- as.data.frame(log10(counts)) %>% 
  utils::stack() %>% 
  mutate(Genotype=samples$Genotype[match(ind,samples$SampleID)],
         TechRep=samples$TechRep[match(ind,samples$SampleID)])

ggplot(dat,aes(x=values,group=ind,col=Genotype)) + 
  geom_density(na.rm = TRUE) + 
  ggtitle("sample raw counts distribution") +
  scale_x_continuous(name="per gene raw counts (log10)") + 
  theme_bw()

ggplot(dat,aes(x=values,group=ind,col=TechRep)) + 
  geom_density(na.rm = TRUE) + 
  ggtitle("sample raw counts distribution") +
  scale_x_continuous(name="per gene raw counts (log10)") + 
  theme_bw()

#' `r emoji("point_right")` **Samples have different sequencing depth that is in agreement with the presence of technical replicates**
#' 
#' * Export raw expression data
dir.create(here("data/analysis/salmon"),showWarnings=FALSE,recursive=TRUE)
write.csv(counts,file=here("data/analysis/salmon/raw-unormalised-gene-expression_data.csv"))
#' 
#' <hr />
#' &nbsp;
#' 
#' # Data normalisation 
#' ## Preparation
#' 
#' For visualization, the data is submitted to a variance stabilization
#' transformation using _DESeq2_. The dispersion is estimated independently
#' of the sample tissue and replicate. 
#'  
dds <- DESeqDataSetFromTximport(
  txi=txi,
  colData = samples,
  design = ~ Genotype)

colnames(dds) <- samples$SampleID

save(dds,file=here("data/analysis/salmon/dds.rda"))

#' ## size factors 
#' (_i.e._ the sequencing library size effect)
#' 
dds <- estimateSizeFactors(dds) %>% 
  suppressMessages()

boxplot(normalizationFactors(dds),
        main="Sequencing libraries size factor",
        las=2,log="y")
abline(h=1, col = "Red", lty = 3)

#' and without outliers:
boxplot(normalizationFactors(dds),
        main="Sequencing libraries size factor",
        las=2,log="y", outline=FALSE)
abline(h=1, col = "Red", lty = 3)

#' `r emoji("point_right")` **Seven samples had been clearly under-sequenced.**
#'
#' Assess whether there might be a difference in library size linked to a
#' given metadata
boxplot(split(t(normalizationFactors(dds)),dds$Genotype),las=2,
        main="Sequencing libraries size factor by Genotype",
        outline=FALSE)

#' `r emoji("point_right")` **The scaling factor distribution is similar for all genotypes with line3 having less lowly sequenced technical replicates.**
#'
plot(colMeans(normalizationFactors(dds)),
     log10(colSums(counts(dds))),ylab="log10 raw depth",
     xlab="average scaling factor",
     col=rainbow(n=nlevels(dds$Genotype))[as.integer(dds$Genotype)],pch=19)
legend("bottomright",fill=rainbow(n=nlevels(dds$Genotype)),
       legend=levels(dds$Genotype),cex=0.6)

#' `r emoji("point_right")` **The scaling factor appears linear for the samples with sufficient sequencing depth. Otherwise, not unintuitively, the distribution is logarithmic**
#'
#' ## Variance Stabilising Transformation
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
vst <- vst - min(vst)

#' ## Validation
#' 
#' let's look at standard deviations before and after VST normalization. 
#' This plot is to see whether there is a dependency of SD on the mean. 
#' 
#' Before:  
meanSdPlot(log1p(counts(dds))[rowSums(counts(dds))>0,])

#' After:
meanSdPlot(vst[rowSums(vst)>0,])

#' After VST normalization, the red line is almost horizontal which indicates
#' no dependency of variance on mean (homoscedastic).
#' 
#' `r emoji("point_right")` **We can conclude that the variance stabilization worked adequately**
#' 
#' <hr />
#' &nbsp;
#' 
#' ## QC on the normalised data
#' ### PCA
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)

#' Using PCAtools
p <- pca(vst,colData(dds))

#' ### Scree plot
#' 
#' We define the number of variable of the model:
nvar <- 1

#' An the number of possible combinations
nlevel <- nlevels(dds$Genotype)

#' We devise the optimal number of components using two methods
horn <- suppressWarnings(parallelPCA(vst))
elbow <- findElbowPoint(p$variance)

#' We plot the percentage explained by different components and try to empirically assess whether
#' the observed number of components would be in agreement with our model's assumptions.
#' 
#' * the red line represents number of variables in the model  
#' * the orange line represents number of variable combinations.
#' * the black dotted, annotate lines represent the optimal number of components 
#' reported by the horn and elbow methods.
#' 
ggplot(tibble(x=1:length(percent),y=cumsum(percent),p=percent),aes(x=x,y=y)) +
  geom_line() + geom_col(aes(x,p)) + scale_y_continuous("variance explained (%)",limits=c(0,100)) +
  scale_x_continuous("Principal component",breaks=1:length(percent),minor_breaks=NULL) + 
  geom_vline(xintercept=nvar,colour="red",linetype="dashed",linewidth=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nvar],colour="red",linetype="dashed",linewidth=0.5) +
  geom_vline(xintercept=nlevel,colour="orange",linetype="dashed",linewidth=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nlevel],colour="orange",linetype="dashed",linewidth=0.5) +
  geom_vline(xintercept=c(horn$n,elbow),colour="black",linetype="dotted",linewidth=0.5) +
  geom_label(aes(x = horn$n + 1, y = cumsum(percent)[horn$n],label = 'Horn', vjust = 1)) +
  geom_label(aes(x = elbow + 1, y = cumsum(percent)[elbow],label = 'Elbow', vjust = 1))

#' `r emoji("point_right")` **The first component explains 50% of the data variance. Both metrics, Horn and Elbow suggest that three or four components are those that are informative. Indeed the slope of the curve is fairly linear past PC3 and that would indicate that the remaining PCs only capture sample specific noise. While this is only empirical, the scree plot support having only few variables of importance in the dataset.**
#'

#' ### PCA plot
pc.dat <- bind_cols(PC1=pc$x[,1],
                    PC2=pc$x[,2],
                    PC3=pc$x[,3],
                    as.data.frame(colData(dds)))

#' * PC1 vs PC2
p1 <- ggplot(pc.dat,aes(x=PC1,y=PC2,col=Genotype,shape=BioRep,text=SampleID)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")

ggplotly(p1) %>% 
  layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC2 (",percent[2],"%)",sep=""))) %>% suppressWarnings()

#' The same as a biplot
biplot(p,lab=samples$BioRep,
       colby = 'Genotype',
       colLegendTitle = 'Genotype',
       encircle = TRUE,
       encircleFill = TRUE,
       legendPosition = 'top', 
       legendLabSize = 16, legendIconSize = 8.0)

#' `r emoji("point_right")` **The first dimension separates the lines from WT (T89). The second dimension separates the lines. The technical replicates cluster close together.**
#'
#' * PC1 vs PC3
p2 <- ggplot(pc.dat,aes(x=PC1,y=PC3,col=Genotype,shape=BioRep,text=SampleID)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")

ggplotly(p2) %>% 
  layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC3 (",percent[3],"%)",sep=""))) %>% suppressWarnings()

#' The same as a biplot
biplot(p,x = 'PC1', y = 'PC3',
       lab = samples$BioRep,
       colby = 'Genotype',
       colLegendTitle = 'Genotype',
       encircle = TRUE,
       encircleFill = TRUE,
       legendPosition = 'top', 
       legendLabSize = 16, legendIconSize = 8.0)

#' `r emoji("point_right")` **Interestingly, the third dimension separates biological replicates from line3. Note that this only explains 8% of the data variance.**
#"
#' ```{r subplot, out.width = '100%'}
#' subplot(style(p1, showlegend = FALSE), p2,
#'         titleX = TRUE, titleY = TRUE, nrows = 1, margin = c(0.05, 0.05, 0, 0))
#' ```

#' ### Pairs plot
#' This allows for looking at more dimensions, five by default
#' 
suppressMessages(pairsplot(p,colby='Genotype',shape='BioRep'))

#' `r emoji("point_right")` **PC1 and PC2 are explained by the genotypes. PC3 to PC5, each explain the variance within one genotype. **
#'
#' ### Loadings
#' Loadings, _i.e._ the individual effect of every gene in a component can be studied. Here the most important ones are visualized throughout the different PCs
plotloadings(p,
             rangeRetain = 0.01,
             labSize = 4.0,
             title = 'Loadings plot',
             subtitle = 'PC1 to PC5',
             caption = 'Top 1% variables',
             drawConnectors = TRUE)

#' `r emoji("point_right")` **Among the nine genes selected that have the most important loading to the first five PCs, none of them is in the list of gene of interest. It might be interesting to look these genes up.**
#' 
#' ### Samples Distance
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(sampleDistMatrix) <- paste(dds$Genotype,dds$BioRep)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=pal)

#' `r emoji("point_right")` **The samples group be genotype. Not all the technical replicates group as expected, but that is not unexpected given the shallow sequencing depth of some of these.**
#'
#' ## Sequencing depth
#' The figures show the number of genes expressed per condition at different expression cutoffs. The scale on the lower plot is the same as on the upper.
#' The first plot is a heatmap showing the number of genes above a given cutoff. The second plot shows it as a ratio of the number of genes expressed for (a)
#' given variable(s) divided by the average number of genes expressed at that cutoff across all variable(s). The latter plot is of course biased at higher cutoff 
#' as the number of genes becomes smaller and smaller.
#' The point of these two plots is to assert whether the number of genes expressed varies between conditions, as this would break some assumptions for normalisation and differential expression.
conds <- dds$Genotype
dev.null <- rangeSamplesSummary(counts=vst,
                                conditions=conds,
                                nrep=3)

#' `r emoji("point_right")` **There is almost no difference in the number of genes expressed by the different genotypes.**
#'
#' Plotting the number of genes that are expressed (at any level)
do.call(rbind,split(t(nrow(vst) - colSums(vst==0)),samples$Genotype)) %>% as.data.frame() %>% 
  rownames_to_column("Genotype") %>% pivot_longer(starts_with("V")) %>% 
  ggplot(aes(x=Genotype, y=value, fill=Genotype)) + geom_dotplot(binaxis = "y", stackdir = "center")

#' `r emoji("point_right")` **Here the results are clearly affected by the lower sequencing of some of the samples.**
#'
#' ## Heatmap
#' 
#' Here we want to visualise all the informative genes as a heatmap. We first filter the genes to remove those below the selected noise/signal cutoff. 
#' The method employed here is naive, and relies on observing a sharp decrease in the number of genes within the first few low level of expression. 
#' Using an independent filtering method, such as implemented in DESeq2 would be more accurate, but for the purpose of QA validation, a naive approach is sufficient.
#' Note that a sweet spot for computation is around 20 thousand genes, as building the hierarchical clustering for the heatmap scales non-linearly.
#'
sels <- rangeFeatureSelect(counts=vst,
                           conditions=conds,
                           nrep=3) %>% 
  suppressWarnings()

#' `r emoji("point_right")` **Here a cutoff of 1 is applied, as a naive way to separate the signal from the noise.**
vst.cutoff <- 1

nn <- nrow(vst[sels[[vst.cutoff+1]],])
tn <- nrow(vst)
pn <- round(nn * 100/ tn, digits=1)

#' `r emoji("warning")` **`r pn`%** (`r nn`) of total `r tn` genes are plotted below:
#' 
mat <- t(scale(t(vst[sels[[vst.cutoff+1]],])))
hm <- pheatmap(mat,
               color = hpal,
               border_color = NA,
               clustering_distance_rows = "correlation",
               clustering_method = "ward.D2",
               show_rownames = FALSE,
               labels_col = conds,
               angle_col = 90,
               legend = FALSE)

#' `r emoji("point_right")` **The heatmap shows the expected clustering.**
#'
#' ## Clustering of samples
#'
#' Below we assess the previous dendrogram's reproducibility and plot the clusters with au and bp where:
#' 
#' * __au (Approximately Unbiased): computed by multiscale bootstrap resampling__ `r emoji("point_left")` a better approximation
#' * __bp (Bootstrap Probability): computed by normal bootstrap resampling__
#' 
#' `r emoji("warning")`Clusters with AU larger than 95% are highlighted by rectangles, which are strongly supported by data
#' 
hm.pvclust <- pvclust(data = t(scale(t(vst[sels[[vst.cutoff+1]],]))),
                       method.hclust = "ward.D2", 
                       nboot = 100, parallel = TRUE)

plot(hm.pvclust, labels = paste(dds$Genotype,dds$BioRep))
pvrect(hm.pvclust)

#' `r emoji("point_right")` **The clustering results shows after bootstrapping that the technical replicates are all as they should, closest with each-other.**
#'
#' <details><summary>bootstrapping results as a table</summary>
#' ```{r bootstrapping results as a table}
#' print(hm.pvclust, digits=3)
#' ```
#' </details>
#' 
#' # Merging technical replicates
samples$BioID <- paste0(samples$Genotype,"_",samples$BioRep)
txi$counts <- sapply(split.data.frame(t(txi$counts),samples$BioID),colSums)
txi$length <- sapply(split.data.frame(t(txi$length),samples$BioID),colMaxs)
samples <- samples[match(colnames(txi$counts),samples$BioID),]
counts <- txi$counts

#' ## Recreate the dds
dds <- DESeqDataSetFromTximport(
  txi=txi,
  colData = samples,
  design = ~ Genotype)

#' ## QA again (quick)
#' * "Not expressed" genes
sel <- rowSums(counts) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(counts),digits=1),
        sum(sel),
        nrow(counts))

#' ### Sequencing depth
#' * Let us take a look at the sequencing depth, colouring by CHANGEME
dat <- tibble(x=colnames(counts),y=colSums(counts)) %>% 
  bind_cols(samples)

ggplot(dat,aes(x,y,fill=Genotype)) + 
  geom_col() + 
  scale_y_continuous(name="reads") +
  facet_grid(~ factor(Genotype), scales = "free") +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=90,size=4),
        axis.title.x=element_blank())

#' `r emoji("point_right")` **The raw sequencing depth is now very similar between replicates.**
#'
#' ### per-gene mean expression
#' 
#' _i.e._ the mean raw count of every gene across samples is calculated
#' and displayed on a log10 scale.
#' 

ggplot(data.frame(value=log10(rowMeans(counts))),aes(x=value)) + 
  geom_density(na.rm = TRUE) +
  ggtitle("gene mean raw counts distribution") +
  scale_x_continuous(name="mean raw counts (log10)") + 
  theme_bw()

#' `r emoji("point_right")` **The cumulative gene coverage is as expected, the great sequencing depth means that we have a broad left shoulder, likely indicative of an elevated signal-to-noise ratio.**
#' 
#' ### Per-sample expression
dat <- as.data.frame(log10(counts)) %>% 
  utils::stack() %>% 
  mutate(Genotype=samples$Genotype[match(ind,samples$BioID)])

ggplot(dat,aes(x=values,group=ind,col=Genotype)) + 
  geom_density(na.rm = TRUE) + 
  ggtitle("sample raw counts distribution") +
  scale_x_continuous(name="per gene raw counts (log10)") + 
  theme_bw()

#' `r emoji("point_right")` **All samples have a similar profile.**
#' 
#' * Export raw expression data
write.csv(counts,file=here("data/analysis/salmon/raw-unormalised_technical-replicates_combined-gene-expression_data.csv"))

#' ### Normalisation
dds <- estimateSizeFactors(dds) %>% 
  suppressMessages()

boxplot(normalizationFactors(dds),
        main="Sequencing libraries size factor",
        las=2,log="y")
abline(h=1, col = "Red", lty = 3)

#' and without outliers:
boxplot(normalizationFactors(dds),
        main="Sequencing libraries size factor",
        las=2,log="y", outline=FALSE)
abline(h=1, col = "Red", lty = 3)

#' `r emoji("point_right")` **The difference in sequencing depth varies between -15/+25%, which is little.**
#'
#' Assess whether there might be a difference in library size linked to a
#' given metadata
boxplot(split(t(normalizationFactors(dds)),dds$Genotype),las=2,
        main="Sequencing libraries size factor by Genotype",
        outline=FALSE)

#' `r emoji("point_right")` **The scaling factor distribution varies between genotypes, possibly indicative that different number of genes are expressed in the different genotypes.**
#'
plot(colMeans(normalizationFactors(dds)),
     log10(colSums(counts(dds))),ylab="log10 raw depth",
     xlab="average scaling factor",
     col=rainbow(n=nlevels(dds$Genotype))[as.integer(dds$Genotype)],pch=19)
legend("bottomright",fill=rainbow(n=nlevels(dds$Genotype)),
       legend=levels(dds$Genotype),cex=0.6)

#' `r emoji("point_right")` **If we assume a linear relationship across the "line" samples, then T89 would seem to have a lower scaling factor for a higher depth, a sign that the number of genes expressed in T89 (WT) would be different from than in the mutated lines.**
#'
#' ### Variance Stabilising Transformation
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
vst <- vst - min(vst)

#' #### Validation
#' 
#' let's look at standard deviations before and after VST normalization. 
#' This plot is to see whether there is a dependency of SD on the mean. 
#' 
#' Before:  
meanSdPlot(log1p(counts(dds))[rowSums(counts(dds))>0,])

#' After:
meanSdPlot(vst[rowSums(vst)>0,])

#' `r emoji("point_right")` **The VST is working as expected**
#' 
#' ### PCA
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)

#' Using PCAtools
p <- pca(vst,colData(dds))

#' #### Scree plot
#' We devise the optimal number of components using two methods
horn <- suppressWarnings(parallelPCA(vst))
elbow <- findElbowPoint(p$variance)

#' We plot the percentage explained by different components and try to empirically assess whether
#' the observed number of components would be in agreement with our model's assumptions.
#' 
#' * the red line represents number of variables in the model  
#' * the orange line represents number of variable combinations.
#' * the black dotted, annotate lines represent the optimal number of components 
#' reported by the horn and elbow methods.
#' 
ggplot(tibble(x=1:length(percent),y=cumsum(percent),p=percent),aes(x=x,y=y)) +
  geom_line() + geom_col(aes(x,p)) + scale_y_continuous("variance explained (%)",limits=c(0,100)) +
  scale_x_continuous("Principal component",breaks=1:length(percent),minor_breaks=NULL) + 
  geom_vline(xintercept=nvar,colour="red",linetype="dashed",linewidth=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nvar],colour="red",linetype="dashed",linewidth=0.5) +
  geom_vline(xintercept=nlevel,colour="orange",linetype="dashed",linewidth=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nlevel],colour="orange",linetype="dashed",linewidth=0.5) +
  geom_vline(xintercept=c(horn$n,elbow),colour="black",linetype="dotted",linewidth=0.5) +
  geom_label(aes(x = horn$n + 1, y = cumsum(percent)[horn$n],label = 'Horn', vjust = 1)) +
  geom_label(aes(x = elbow + 1, y = cumsum(percent)[elbow],label = 'Elbow', vjust = 1))

#' `r emoji("point_right")` **The first component explains 53% of the data variance. Both metrics, Horn and Elbow suggest that two or three components are those that are informative. Indeed the slope of the curve is fairly linear past PC3 and that would indicate that the remaining PCs only capture sample specific noise. While this is only empirical, the scree plot support having only few variables of importance in the dataset.**
#'

#' #### PCA plot
pc.dat <- bind_cols(PC1=pc$x[,1],
                    PC2=pc$x[,2],
                    PC3=pc$x[,3],
                    as.data.frame(colData(dds)))

#' * PC1 vs PC2
p1 <- ggplot(pc.dat,aes(x=PC1,y=PC2,col=Genotype,shape=BioRep,text=SampleID)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")

ggplotly(p1) %>% 
  layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC2 (",percent[2],"%)",sep=""))) %>% suppressWarnings()

#' The same as a biplot
biplot(p,lab=samples$BioRep,
       colby = 'Genotype',
       colLegendTitle = 'Genotype',
       encircle = TRUE,
       encircleFill = TRUE,
       legendPosition = 'top', 
       legendLabSize = 16, legendIconSize = 8.0)

#' `r emoji("point_right")` **The first dimension separates the lines from WT (T89). The second dimension separates the lines. The technical replicates cluster close together.**
#'
#' * PC1 vs PC3
p2 <- ggplot(pc.dat,aes(x=PC1,y=PC3,col=Genotype,shape=BioRep,text=SampleID)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")

ggplotly(p2) %>% 
  layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC3 (",percent[3],"%)",sep=""))) %>% suppressWarnings()

#' The same as a biplot
biplot(p,x = 'PC1', y = 'PC3',
       lab = samples$BioRep,
       colby = 'Genotype',
       colLegendTitle = 'Genotype',
       encircle = TRUE,
       encircleFill = TRUE,
       legendPosition = 'top', 
       legendLabSize = 16, legendIconSize = 8.0)

#' `r emoji("point_right")` **Interestingly, the third dimension separates biological replicates from line3. Note that this only explains 8.52% of the data variance.**
#"
#' ```{r subplot2, out.width = '100%'}
#' subplot(style(p1, showlegend = FALSE), p2,
#'         titleX = TRUE, titleY = TRUE, nrows = 1, margin = c(0.05, 0.05, 0, 0))
#' ```

#' ### Pairs plot
#' This allows for looking at more dimensions, five by default
#' 
suppressMessages(pairsplot(p,colby='Genotype',shape='BioRep'))

#' `r emoji("point_right")` **PC1 and PC2 are explained by the genotypes. PC3 to PC4, each explain the variance within one "line" genotype. PC5 is harder to explain, but seem influenced by some replicates.**
#'
#' ### Loadings
#' Loadings, _i.e._ the individual effect of every gene in a component can be studied. Here the most important ones are visualized throughout the different PCs
plotloadings(p,
             rangeRetain = 0.01,
             labSize = 4.0,
             title = 'Loadings plot',
             subtitle = 'PC1 to PC5',
             caption = 'Top 1% variables',
             drawConnectors = TRUE)

#' `r emoji("point_right")` **Among the eight genes selected that have the most important loading to the first five PCs, none of them is in the list of gene of interest. It might be interesting to look these genes up.**
#' 
#' ### Samples Distance
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(sampleDistMatrix) <- paste(dds$Genotype,dds$BioRep)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=pal)

#' `r emoji("point_right")` **The samples group be genotype.**
#'
#' ## Sequencing depth
#' The figures show the number of genes expressed per condition at different expression cutoffs. The scale on the lower plot is the same as on the upper.
#' The first plot is a heatmap showing the number of genes above a given cutoff. The second plot shows it as a ratio of the number of genes expressed for (a)
#' given variable(s) divided by the average number of genes expressed at that cutoff across all variable(s). The latter plot is of course biased at higher cutoff 
#' as the number of genes becomes smaller and smaller.
#' The point of these two plots is to assert whether the number of genes expressed varies between conditions, as this would break some assumptions for normalisation and differential expression.
conds <- dds$Genotype
dev.null <- rangeSamplesSummary(counts=vst,
                                conditions=conds,
                                nrep=3)

#' `r emoji("point_right")` **There is almost no difference in the number of genes expressed by the different genotypes at lower level of expression, but T89 has more genes expressed and at a highly level of expression.**
#'
#' Plotting the number of genes that are expressed (at any level)
do.call(rbind,split(t(nrow(vst) - colSums(vst==0)),samples$Genotype)) %>% as.data.frame() %>% 
  rownames_to_column("Genotype") %>% pivot_longer(starts_with("V")) %>% 
  ggplot(aes(x=Genotype, y=value, fill=Genotype)) + geom_dotplot(binaxis = "y", stackdir = "center")

#' `r emoji("point_right")` **This is interesting, overall T89 has less genes expressed but these are expressed at a higher level.**
#'
#' ## Heatmap
#' 
#' Here we want to visualise all the informative genes as a heatmap. We first filter the genes to remove those below the selected noise/signal cutoff. 
#' The method employed here is naive, and relies on observing a sharp decrease in the number of genes within the first few low level of expression. 
#' Using an independent filtering method, such as implemented in DESeq2 would be more accurate, but for the purpose of QA validation, a naive approach is sufficient.
#' Note that a sweet spot for computation is around 20 thousand genes, as building the hierarchical clustering for the heatmap scales non-linearly.
#'
sels <- rangeFeatureSelect(counts=vst,
                           conditions=conds,
                           nrep=3) %>% 
  suppressWarnings()

#' `r emoji("point_right")` **Here a cutoff of 1 is applied, as a naive way to separate the signal from the noise.**
vst.cutoff <- 1

nn <- nrow(vst[sels[[vst.cutoff+1]],])
tn <- nrow(vst)
pn <- round(nn * 100/ tn, digits=1)

#' `r emoji("warning")` **`r pn`%** (`r nn`) of total `r tn` genes are plotted below:
#' 
mat <- t(scale(t(vst[sels[[vst.cutoff+1]],])))
hm <- pheatmap(mat,
               color = hpal,
               border_color = NA,
               clustering_distance_rows = "correlation",
               clustering_method = "ward.D2",
               show_rownames = FALSE,
               labels_col = conds,
               angle_col = 90,
               legend = FALSE)

#' `r emoji("point_right")` **The heatmap shows the expected clustering and a lot of potential for differential expression.**
#'
#' ## Clustering of samples
#'
#' Below we assess the previous dendrogram's reproducibility and plot the clusters with au and bp where:
#' 
#' * __au (Approximately Unbiased): computed by multiscale bootstrap resampling__ `r emoji("point_left")` a better approximation
#' * __bp (Bootstrap Probability): computed by normal bootstrap resampling__
#' 
#' `r emoji("warning")`Clusters with AU larger than 95% are highlighted by rectangles, which are strongly supported by data
#' 
hm.pvclust <- pvclust(data = t(scale(t(vst[sels[[vst.cutoff+1]],]))),
                      method.hclust = "ward.D2", 
                      nboot = 100, parallel = TRUE)

plot(hm.pvclust, labels = paste(dds$Genotype,dds$BioRep))
pvrect(hm.pvclust)

#'
#' <hr />
#' &nbsp;
#' 
#' # Summary
#' `r emoji("star")` **The data is of good quality. The initial data had multiple technical replicates that were merged after an initial assessment.**
#' 
#' `r emoji("star")` **The technical replicate merged data is also of good quality, with all assessment plots separating the genotypes and grouping the replicates as they should.**
#' 
#' `r emoji("star")` **The wild-type T89 seem to have a less diversified transcriptome - less genes expressed - at a comparatively higher level. This might be an artifact of the normalization that expects the same number of genes to be expressed across all genotypes.**
#' 
#' ```{r empty,eval=FALSE,echo=FALSE}
#' ```
#'
#' # Session Info
#' <details><summary>Session Info</summary>
#' ```{r session info}
#' sessionInfo()
#' ```
#' </details>
#'   
#' &nbsp;
#' 
