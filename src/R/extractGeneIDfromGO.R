### Exportation of DEGs requested by Daniela
library(DESeq2)
load(here("analysis/salmon/dds_mergeTechRep.rda"))
dds$Genotype <- relevel(dds$Genotype,ref = "T89")

vsd <- varianceStabilizingTransformation(dds,blind=FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)

source(here("UPSCb-common/src/R/featureSelection.R"))
source(here("UPSCb-common/src/R/topGoUtilities.R"))

background <- rownames(vst)[featureSelect(vst,dds$Genotype,exp=0.00001)]
goannot <- prepAnnot(mapping = "/mnt/picea/storage/reference/Populus-tremula/v2.2/gopher/gene_to_go.tsv")
PvalCutoff = 0.01
PadjustMethod = "BH"

res.list <- list("CommonLine26andLine03"=list("all"=read.csv(here("analysis/DE/commonL26L3_all.txt"),header = F)$V1,
																							"up"=read.csv(here("analysis/DE/commonL26L3_up.txt"),header = F)$V1,
																							"dn"=read.csv(here("analysis/DE/commonL26L3_down.txt"),header = F)$V1),
								 "ExclusivelyLine26"=list("all"=read.csv(here("analysis/DE/specificL26_all.txt"),header = F)$V1,
								 												 "up"=read.csv(here("analysis/DE/specificL26_up.txt"),header = F)$V1,
								 												 "dn"=read.csv(here("analysis/DE/specificL26_down.txt"),header = F)$V1),
								 "ExclusivelyLine03"=list("all"=read.csv(here("analysis/DE/specificL3_all.txt"),header = F)$V1,
								 												 "up"=read.csv(here("analysis/DE/specificL3_up.txt"),header = F)$V1,
								 												 "dn"=read.csv(here("analysis/DE/specificL3_down.txt"),header = F)$V1))

enr.list <- lapply(res.list,function(r){
	lapply(r,topGO,background=background,annotation=goannot,alpha=PvalCutoff,p.adjust=PadjustMethod,ontology="BP",getgenes=TRUE)
})

goannottable <- read.table("/mnt/picea/storage/reference/Populus-tremula/v2.2/gopher/gene_to_go.tsv")
colnames(goannottable) <- c("GeneID","GO")

library(xlsx)
#should have do this as lapply

datg <- data.frame(res.list$CommonLine26andLine03$all)
dat <- data.frame(enr.list$CommonLine26andLine03$all$BP)
dat$allgenes <- NULL
for(i in 1:length(dat$siggenes)){
	dat$siggenes[i] <- paste(unlist(dat$siggenes[i]),collapse = "|")
}
write.xlsx(dat, file=here("analysis/topGO/CommonLine26andLine03.xlsx"), sheetName="all_BP", row.names=FALSE)
datg <- data.frame(res.list$CommonLine26andLine03$up)
dat <- data.frame(enr.list$CommonLine26andLine03$up$BP)
dat$allgenes <- NULL
for(i in 1:length(dat$siggenes)){
	dat$siggenes[i] <- paste(unlist(dat$siggenes[i]),collapse = "|")
}
write.xlsx(dat, file=here("analysis/topGO/CommonLine26andLine03.xlsx"), sheetName="up_BP", append=TRUE, row.names=FALSE)
datg <- data.frame(res.list$CommonLine26andLine03$dn)
dat <- data.frame(enr.list$CommonLine26andLine03$dn$BP)
dat$allgenes <- NULL
for(i in 1:length(dat$siggenes)){
	dat$siggenes[i] <- paste(unlist(dat$siggenes[i]),collapse = "|")
}
write.xlsx(dat, file=here("analysis/topGO/CommonLine26andLine03.xlsx"), sheetName="dn_BP", append=TRUE, row.names=FALSE)
dat <- data.frame(res.list$CommonLine26andLine03$all)
colnames(dat) <- "GeneID"
datgo <- left_join(dat,goannottable)
write.xlsx(datgo, file=here("analysis/topGO/CommonLine26andLine03.xlsx"), sheetName="all_gene", append=TRUE, row.names=FALSE)
dat <- data.frame(res.list$CommonLine26andLine03$up)
colnames(dat) <- "GeneID"
datgo <- left_join(dat,goannottable)
write.xlsx(datgo, file=here("analysis/topGO/CommonLine26andLine03.xlsx"), sheetName="up_gene", append=TRUE, row.names=FALSE)
dat <- data.frame(res.list$CommonLine26andLine03$dn)
colnames(dat) <- "GeneID"
datgo <- left_join(dat,goannottable)
write.xlsx(datgo, file=here("analysis/topGO/CommonLine26andLine03.xlsx"), sheetName="dn_gene", append=TRUE, row.names=FALSE)



datg <- data.frame(res.list$ExclusivelyLine26$all)
dat <- data.frame(enr.list$ExclusivelyLine26$all$BP)
dat$allgenes <- NULL
for(i in 1:length(dat$siggenes)){
	dat$siggenes[i] <- paste(unlist(dat$siggenes[i]),collapse = "|")
}
write.xlsx(dat, file=here("analysis/topGO/ExclusivelyLine26.xlsx"), sheetName="all_BP", row.names=FALSE)
# datg <- data.frame(res.list$ExclusivelyLine26$up)
# dat <- data.frame(enr.list$ExclusivelyLine26$up$BP)
# dat$allgenes <- NULL
# for(i in 1:length(dat$siggenes)){
# 	dat$siggenes[i] <- paste(unlist(dat$siggenes[i]),collapse = "|")
# }
# write.xlsx(dat, file=here("analysis/topGO/ExclusivelyLine26.xlsx"), sheetName="up_BP", append=TRUE, row.names=FALSE)
datg <- data.frame(res.list$ExclusivelyLine26$dn)
dat <- data.frame(enr.list$ExclusivelyLine26$dn$BP)
dat$allgenes <- NULL
for(i in 1:length(dat$siggenes)){
	dat$siggenes[i] <- paste(unlist(dat$siggenes[i]),collapse = "|")
}
write.xlsx(dat, file=here("analysis/topGO/ExclusivelyLine26.xlsx"), sheetName="dn_BP", append=TRUE, row.names=FALSE)
dat <- data.frame(res.list$ExclusivelyLine26$all)
colnames(dat) <- "GeneID"
datgo <- left_join(dat,goannottable)
write.xlsx(datgo, file=here("analysis/topGO/ExclusivelyLine26.xlsx"), sheetName="all_gene", append=TRUE, row.names=FALSE)
dat <- data.frame(res.list$ExclusivelyLine26$up)
colnames(dat) <- "GeneID"
datgo <- left_join(dat,goannottable)
write.xlsx(datgo, file=here("analysis/topGO/ExclusivelyLine26.xlsx"), sheetName="up_gene", append=TRUE, row.names=FALSE)
dat <- data.frame(res.list$ExclusivelyLine26$dn)
colnames(dat) <- "GeneID"
datgo <- left_join(dat,goannottable)
write.xlsx(datgo, file=here("analysis/topGO/ExclusivelyLine26.xlsx"), sheetName="dn_gene", append=TRUE, row.names=FALSE)


datg <- data.frame(res.list$ExclusivelyLine03$all)
dat <- data.frame(enr.list$ExclusivelyLine03$all$BP)
dat$allgenes <- NULL
for(i in 1:length(dat$siggenes)){
	dat$siggenes[i] <- paste(unlist(dat$siggenes[i]),collapse = "|")
}
write.xlsx(dat, file=here("analysis/topGO/ExclusivelyLine03.xlsx"), sheetName="all_BP", row.names=FALSE)
datg <- data.frame(res.list$ExclusivelyLine03$up)
dat <- data.frame(enr.list$ExclusivelyLine03$up$BP)
dat$allgenes <- NULL
for(i in 1:length(dat$siggenes)){
	dat$siggenes[i] <- paste(unlist(dat$siggenes[i]),collapse = "|")
}
write.xlsx(dat, file=here("analysis/topGO/ExclusivelyLine03.xlsx"), sheetName="up_BP", append=TRUE, row.names=FALSE)
datg <- data.frame(res.list$ExclusivelyLine03$dn)
dat <- data.frame(enr.list$ExclusivelyLine03$dn$BP)
dat$allgenes <- NULL
for(i in 1:length(dat$siggenes)){
	dat$siggenes[i] <- paste(unlist(dat$siggenes[i]),collapse = "|")
}
write.xlsx(dat, file=here("analysis/topGO/ExclusivelyLine03.xlsx"), sheetName="dn_BP", append=TRUE, row.names=FALSE)
dat <- data.frame(res.list$ExclusivelyLine03$all)
colnames(dat) <- "GeneID"
datgo <- left_join(dat,goannottable)
write.xlsx(datgo, file=here("analysis/topGO/ExclusivelyLine03.xlsx"), sheetName="all_gene", append=TRUE, row.names=FALSE)
dat <- data.frame(res.list$ExclusivelyLine03$up)
colnames(dat) <- "GeneID"
datgo <- left_join(dat,goannottable)
write.xlsx(datgo, file=here("analysis/topGO/ExclusivelyLine03.xlsx"), sheetName="up_gene", append=TRUE, row.names=FALSE)
dat <- data.frame(res.list$ExclusivelyLine03$dn)
colnames(dat) <- "GeneID"
datgo <- left_join(dat,goannottable)
write.xlsx(datgo, file=here("analysis/topGO/ExclusivelyLine03.xlsx"), sheetName="dn_gene", append=TRUE, row.names=FALSE)


