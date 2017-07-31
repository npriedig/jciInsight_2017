################################################################################################
# set working directory to userDirectory and run script below in Rstudio
# ensure required packages are installed (i.e. run installRequiredPackages.R)
# note: scripts involving figures w/ controlled access data (i.e. METABRIC, variants) 
# are not represented here, but should be reproducible w/ methods in manuscript. 
################################################################################################

# working directory MUST be ecRNAseq_brcaBoneMets folder before running script
userDirectory <- "~/Desktop/ecRNAseq_brcaBoneMets"
setwd(userDirectory)
requiredPackages.cran <- c(
  "devtools", 
  "reshape",  
  "RColorBrewer", 
  "ggplot2", 
  "gplots",
  "corrplot",
  "lattice",
  "survminer")
requiredPackages.bioconductor <- c(
  "genefu",
  "ComplexHeatmap",
  "edgeR",
  "biomaRt",
  "DESeq2")
lapply(c(requiredPackages.bioconductor,requiredPackages.cran, "QoRTs"), require, character.only = TRUE)
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
################################################################################################
# functions
# http://stackoverflow.com/questions/15253954/replace-multiple-arguments-with-gsub
mgsub <- function(pattern, replacement, x, ...) {
  if (length(pattern)!=length(replacement)) {
    stop("pattern and replacement do not have the same length.")
  }
  result <- x
  for (i in 1:length(pattern)) {
    result <- gsub(pattern[i], replacement[i], result, ...)
  }
  result
}
# https://stackoverflow.com/questions/2547402/is-there-a-built-in-function-for-finding-the-mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
        ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
      print(plots[[1]])
  } else {
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
      for (i in 1:numPlots) {
          matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
          print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
            layout.pos.col = matchidx$col))
      }
    }
}
################################################################################################
# import raw data
load("dataInput/txi.salmon.boneMets.decalFrozenQC.Rda")
load("dataInput/tumorMatch.scores.df.Rda")
load("dataInput/txi.salmon.boneMets.final.Rda")
load("dataInput/sampleTable.Rda")
load("dataInput/bmfs.Rda")
load("dataInput/salmon.counts.balancedPrimaries.pam50.Rda")
load("dataInput/salmon.tpm.balancedPrimaries.pam50.Rda")
load("dataInput/sampleTable.balancedPrimaries.pam50.Rda")
load("dataInput/boneMetSurvival.GSE12276.Rda")
load("dataInput/gcrmaMedians.GSE12276.Rda")
pam50.gl <- scan("dataInput/ensembl.pam50.gl", what = "character")
clinicallyActionable.gl <- scan("dataInput/dgidb_clinicallyActionable.txt", what = "character")
mart<- useDataset("hsapiens_gene_ensembl", 
  useMart("ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org"))
################################################################################################
# FIGURE 1A, Supplemental Table S1
dataDirectory <- "dataInput/qorts-files/"
decoderFile <- "dataInput/qorts-files/qorts-decoder-frozenVsFFPE.txt"
qorts.res <- read.qc.results.data(dataDirectory,
            decoder.files = decoderFile,
            calc.DESeq2 = TRUE, calc.edgeR = TRUE)
outDirectory <- paste("dataOutput/Figure_1/qorts-summaryPDFs-frozenVsFFPE/", sep = "")
dir.create(outDirectory)
makeMultiPlot.colorByGroup(qorts.res,
                  outfile.dir = outDirectory,
                  plot.device.name = "pdf")
get.summary.table(qorts.res, outfile = paste(outDirectory, "S1.table.txt", sep = "/"))
################################################################################################
# FIGURE 1B
frozenVsFFPE.samples <- c("BP29", "BPF29", "BP51", "BPF51", "BP68", "BPF68", "BP72", "BPF72")
salmon.tpm.selected <- txi.salmon.decalFrozenQC$abundance[,frozenVsFFPE.samples]
number.of.samples.tpmThreshold <- apply(salmon.tpm.selected, MARG = 1, 
  FUN = function(x) length(which(x > 0.5)))
keep.index <- as.vector(which(number.of.samples.tpmThreshold > 1))

cts <- txi.salmon.decalFrozenQC$counts[keep.index,frozenVsFFPE.samples]
y <- DGEList(cts)
y <- calcNormFactors(y)
cpm.norm <- cpm(y, normalized.lib.sizes = TRUE)
cpm.norm.final.log2 <- log2(cpm.norm+1)

pdf(file = "dataOutput/Figure_1/frozenVsFFPE-tmm-cpm_correlations.pdf", height = 8, width = 8)
  color.transparent <- adjustcolor("darkblue", alpha.f = 0.2) 
  #########
  case29.frozenVsFFPE <- cpm.norm.final.log2[,c("BP29", "BPF29")]
  FFPE.frozen.cor <- cor.test(case29.frozenVsFFPE[,1],case29.frozenVsFFPE[,2])
  plot(case29.frozenVsFFPE[,1],case29.frozenVsFFPE[,2],pch = 16, col = color.transparent, 
       cex.main = 1.75, font.lab = 2, main = "BP29-Frozen vs BP29-FFPE",xlab = "BP29-FFPE log2normCPM", 
       ylab = "BP29-Frozen log2normCPM", cex.axis = 1.25, at = c(0,5,10,15), 
       xlim = c(0,16), ylim = c(0,16), cex.lab = 1.3)
  abline(0,1,col = "red3", lwd = 2)
  legend("topleft", bty="n", legend=paste("Pearson r =", 
    signif(as.numeric(FFPE.frozen.cor$estimate), digits=3)),text.font = 2, cex = 1.5)
  #########
  case51.frozenVsFFPE <- cpm.norm.final.log2[,c("BP51", "BPF51")]
  FFPE.frozen.cor <- cor.test(case51.frozenVsFFPE[,1],case51.frozenVsFFPE[,2])
  plot(case51.frozenVsFFPE[,1],case51.frozenVsFFPE[,2],pch = 16, col = color.transparent, 
       cex.main = 1.75, font.lab = 2, main = "BP51-Frozen vs BP51-FFPE",xlab = "BP51-FFPE log2normCPM", 
       ylab = "BP51-Frozen log2normCPM", cex.axis = 1.25, at = c(0,5,10,15), 
       xlim = c(0,16), ylim = c(0,16), cex.lab = 1.3)
  abline(0,1,col = "red3", lwd = 2)
  legend("topleft", bty="n", legend=paste("Pearson r =", 
    signif(as.numeric(FFPE.frozen.cor$estimate), digits=3)),text.font = 2, cex = 1.5)
  #########
  case68.frozenVsFFPE <- cpm.norm.final.log2[,c("BP68", "BPF68")]
  FFPE.frozen.cor <- cor.test(case68.frozenVsFFPE[,1],case68.frozenVsFFPE[,2])
  plot(case68.frozenVsFFPE[,1],case68.frozenVsFFPE[,2],pch = 16, col = color.transparent, 
       cex.main = 1.75, font.lab = 2, main = "BP68-Frozen vs BP68-FFPE",xlab = "BP68-FFPE log2normCPM", 
       ylab = "BP68-Frozen log2normCPM", cex.axis = 1.25, at = c(0,5,10,15), 
       xlim = c(0,16), ylim = c(0,16), cex.lab = 1.3)
  abline(0,1,col = "red3", lwd = 2)
  legend("topleft", bty="n", legend=paste("Pearson r =", 
    signif(as.numeric(FFPE.frozen.cor$estimate), digits=3)),text.font = 2, cex = 1.5)
  #########
  case72.frozenVsFFPE <- cpm.norm.final.log2[,c("BP72", "BPF72")]
  FFPE.frozen.cor <- cor.test(case72.frozenVsFFPE[,1],case72.frozenVsFFPE[,2])
  plot(case72.frozenVsFFPE[,1],case72.frozenVsFFPE[,2],pch = 16, col = color.transparent, 
       cex.main = 1.75, font.lab = 2, main = "BP72-Frozen vs BP72-FFPE",xlab = "BP72-FFPE log2normCPM", 
       ylab = "BP72-Frozen log2normCPM", cex.axis = 1.25, at = c(0,5,10,15), 
       xlim = c(0,16), ylim = c(0,16), cex.lab = 1.3)
  abline(0,1,col = "red3", lwd = 2)
  legend("topleft", bty="n", legend=paste("Pearson r =", 
    signif(as.numeric(FFPE.frozen.cor$estimate), digits=3)),text.font = 2, cex = 1.5)
dev.off()
# Supplemental Figure 1 (1 of 2)
cor.matrix <- cor(cpm.norm.final.log2)
pdf(file = "dataOutput/SFigure_1/frozenVsFFPE_corrMatrix.pdf", height = 8, width = 8)
corrplot(cor.matrix , method = "circle", is.corr = FALSE,tl.col = "black",
  tl.cex = 0.8, order = "hclust")
title("Frozen vs. FFPE Correlation Matrix")
dev.off()
################################################################################################
# FIGURE 1C, Supplemental Table S2
dataDirectory <- "dataInput/qorts-files/"
decoderFile <- "dataInput/qorts-files/boneMets-qortsDecoder-decal-vs-nonDecal.txt"

qorts.res <- read.qc.results.data(dataDirectory, decoder.files = decoderFile,
  calc.DESeq2 = TRUE, calc.edgeR = TRUE)

outDirectory <- paste("dataOutput/Figure_1/qorts-summaryPDFs-decalVsnonDecal/", sep = "")
dir.create(outDirectory)

makeMultiPlot.colorByGroup(qorts.res, outfile.dir = outDirectory, plot.device.name = "pdf")
get.summary.table(qorts.res, outfile = "dataOutput/STables/S2.table.txt")
################################################################################################
# FIGURE 1D
nonDecalvsDecal.samples <- c("BM03", "BM03-DECAL", "BM08", "BM08-DECAL", "22M", "22M-DECAL") 
salmon.tpm.selected <- txi.salmon.decalFrozenQC$abundance[,nonDecalvsDecal.samples]
number.of.samples.tpmThreshold <- apply(salmon.tpm.selected, MARG = 1, 
  FUN = function(x) length(which(x > 0.5)))
keep.index <- as.vector(which(number.of.samples.tpmThreshold > 1))

cts <- txi.salmon.decalFrozenQC$counts[keep.index,nonDecalvsDecal.samples]
y <- DGEList(cts)
y <- calcNormFactors(y)
cpm.norm <- cpm(y, normalized.lib.sizes = TRUE)
cpm.norm.final.log2 <- log2(cpm.norm+1)

pdf(file = "dataOutput/Figure_1/nonDecalvsDecal-tmm-cpm_correlations.pdf", height = 8, width = 8) 
  #########
  BM03.frozenVsFFPE <- cpm.norm.final.log2[,c("BM03-DECAL", "BM03")]
  FFPE.frozen.cor <- cor.test(BM03.frozenVsFFPE[,1],BM03.frozenVsFFPE[,2])
  plot(BM03.frozenVsFFPE[,1],BM03.frozenVsFFPE[,2],pch = 16, col = color.transparent, 
    cex.main = 1.75, font.lab = 2, main = "BM03 Nondecalcified vs BM03 Decalcified",
    xlab = "BM03 Decalcified log2normCPM", ylab = "BM03 Nondecalcified log2normCPM",
    cex.axis = 1.25, at = c(0,5,10,15), xlim = c(0,16), ylim = c(0,16), cex.lab = 1.3)
  abline(0,1,col = "red3", lwd = 2)
  legend("topleft", bty="n", legend=paste("Pearson r =", 
    signif(as.numeric(FFPE.frozen.cor$estimate), digits=3)),text.font = 2, cex = 1.5)
  #########
  BM08.frozenVsFFPE <- cpm.norm.final.log2[,c("BM08-DECAL", "BM08")]
  FFPE.frozen.cor <- cor.test(BM08.frozenVsFFPE[,1],BM08.frozenVsFFPE[,2])
  plot(BM08.frozenVsFFPE[,1],BM08.frozenVsFFPE[,2],pch = 16, col = color.transparent, 
    cex.main = 1.75, font.lab = 2, main = "BM08 Nondecalcified vs BM08 Decalcified",
    xlab = "BM08 Decalcified log2normCPM", ylab = "BM08 Nondecalcified log2normCPM",
    cex.axis = 1.25, at = c(0,5,10,15), xlim = c(0,16), ylim = c(0,16), cex.lab = 1.3)
  abline(0,1,col = "red3", lwd = 2)
  legend("topleft", bty="n", legend=paste("Pearson r =", 
    signif(as.numeric(FFPE.frozen.cor$estimate), digits=3)),text.font = 2, cex = 1.5)
  #########
  BM22.frozenVsFFPE <- cpm.norm.final.log2[,c("22M-DECAL", "22M")]
  FFPE.frozen.cor <- cor.test(BM22.frozenVsFFPE[,1],BM22.frozenVsFFPE[,2])
  plot(BM22.frozenVsFFPE[,1],BM22.frozenVsFFPE[,2],pch = 16, col = color.transparent, 
    cex.main = 1.75, font.lab = 2, main = "BM22 Nondecalcified vs BM22 Decalcified",
    xlab = "BM22 Decalcified log2normCPM", ylab = "BM22 Nondecalcified log2normCPM", 
    cex.axis = 1.25, at = c(0,5,10,15), xlim = c(0,16), ylim = c(0,16), cex.lab = 1.3)
  abline(0,1,col = "red3", lwd = 2)
  legend("topleft", bty="n", legend=paste("Pearson r =", 
    signif(as.numeric(FFPE.frozen.cor$estimate), digits=3)),text.font = 2, cex = 1.5)
dev.off()

# Supplemental Figure 1 (2 of 2)
cor.matrix <- cor(cpm.norm.final.log2)
pdf(file = "dataOutput/SFigure_1/nonDecalvsDecal_corrMatrix.pdf", height = 8, width = 8)
corrplot(cor.matrix , method = "circle", is.corr = FALSE,tl.col = "black",tl.cex = 0.8, order = "hclust")
title("Nondecalcified vs. Decalcified Correlation Matrix")
dev.off()
########################################################################################################
# Supplemental Figure 2
dataDirectory <- "dataInput/qorts-files/"
decoderFile <- "dataInput/qorts-files/boneMets-qortsDecoder-finalCohort.txt"
qorts.res <- read.qc.results.data(dataDirectory, decoder.files = decoderFile,
  calc.DESeq2 = TRUE, calc.edgeR = TRUE)
outDirectory <- paste("dataOutput/SFigure_2/qorts-summaryPDFs-finalCohort/", sep = "")
dir.create(outDirectory)
makeMultiPlot.colorBySample(qorts.res, outfile.dir = outDirectory, plot.device.name = "pdf")
get.summary.table(qorts.res, outfile = "dataOutput/STables/S5.table.txt")
########################################################################################################
# Supplemental Figure 3, # tumorMatch function and example workflow in 'additional-scripts'
pdf(file = "dataOutput/SFigure_3/bone-met-tumorMatch.pdf", height = 8, width = 8)
corrplot(as.matrix(tumorMatch.scores.df), method = "square", 
  is.corr = FALSE,tl.col = "black",tl.cex = 0.8)
dev.off()
########################################################################################################
# Figure 2A
salmon.tpm.selected <- txi.salmon$abundance[,rownames(sampleTable)]
number.of.samples.tpmThreshold <- apply(salmon.tpm.selected, MARG = 1, 
  FUN = function(x) length(which(x > 0.5)))
keep.index <- as.vector(which(number.of.samples.tpmThreshold > 3))

salmon.counts.selected <- txi.salmon$counts[keep.index,rownames(sampleTable)]
y <- DGEList(salmon.counts.selected)
y <- calcNormFactors(y)
cpm.norm <- cpm(y, normalized.lib.sizes = TRUE)
cpm.norm.final.log2 <- log2(cpm.norm+1)

decalcified.status <- mgsub(c("Pos","Neg"), c("darkgreen","black"), sampleTable$Decalcified)
PrimaryMet.colors <- mgsub(c("primary","metastasis"), c("blue2", "firebrick3"), sampleTable$Tumor)
clab <- cbind(decalcified.status, PrimaryMet.colors)
colnames(clab)=c("Decalcification", "Tumor")
breaks=seq(-2, 2, by=0.1)
mycol <- colorpanel(n=length(breaks)-1,low="blue3",mid="white",high="red3")
mydist <- function(x) {as.dist(1 - cor(t(x), use = "pa"))}
myclust <- function(x) {hclust(x,method="average")}
geneNumber = round(nrow(cpm.norm.final.log2)/20)
topVarList <- head(order(-apply(cpm.norm.final.log2, FUN = IQR , MARGIN = 1)),geneNumber)
topVarGenes <- row.names(cpm.norm.final.log2[topVarList,])
data.hm <- as.matrix(cpm.norm.final.log2[topVarGenes,])
pdf(file = "dataOutput/Figure_2/boneMets-heatmap.pdf", width=6, height=10)
  heatmap.3(data.hm, 
    distfun = mydist, hclustfun = myclust, 
    scale = "row", dendrogram = "column", ColSideColors=clab,
    col=mycol, cexRow = 0.5, cexCol = 1.0, trace="none", key=FALSE,
    margins=c(5,5), ColSideColorsSize=3,
    keysize=1, key.title = NA, offsetRow = 0, offsetCol = 0, 
    adjCol = 1, breaks = breaks, symkey = TRUE, labRow = "")
dev.off()

########################################################################################################
# SFigure X
sampleTable.mets <- subset(sampleTable, Tumor == "metastasis")
salmon.tpm.selected.mets <- txi.salmon$abundance[,rownames(sampleTable.mets)]
number.of.samples.tpmThreshold.mets <- apply(salmon.tpm.selected.mets, MARG = 1, 
  FUN = function(x) length(which(x > 0.5)))
keep.index.mets <- as.vector(which(number.of.samples.tpmThreshold.mets > 2))

salmon.counts.selected.mets <- txi.salmon$counts[keep.index.mets,
rownames(sampleTable.mets)]
y.mets <- DGEList(salmon.counts.selected.mets)
y.mets <- calcNormFactors(y.mets)
cpm.norm.mets <- cpm(y.mets, normalized.lib.sizes = TRUE)
cpm.norm.final.log2.mets <- log2(cpm.norm.mets+1)

decalcified.status.mets <- mgsub(c("Pos","Neg"), c("darkgreen","black"), 
  sampleTable.mets$Decalcified)
PrimaryMet.colors.mets <- mgsub(c("primary","metastasis"), c("blue2", "firebrick3"), 
  sampleTable.mets$Tumor)
clab.mets <- cbind(decalcified.status.mets, PrimaryMet.colors.mets)
colnames(clab.mets)=c("Decalcification", "Tumor")
geneNumber.mets = round(nrow(cpm.norm.final.log2.mets)/20)
topVarList.mets <- head(order(-apply(cpm.norm.final.log2.mets, 
  FUN = IQR , MARGIN = 1)),geneNumber.mets)
topVarGenes.mets <- row.names(cpm.norm.final.log2.mets[topVarList.mets,])
data.hm.mets <- as.matrix(cpm.norm.final.log2.mets[topVarGenes.mets,])
pdf(file = paste("dataOutput/SFigure_X/metsOnly-heatmap.pdf"), width=4, height=8)
  heatmap.3(data.hm.mets, 
    distfun = mydist, hclustfun = myclust, 
    scale = "row", dendrogram = "column", ColSideColors=clab.mets,
    col=mycol, cexRow = 0.5, cexCol = 1.0, trace="none", key=FALSE,
    margins=c(5,5), ColSideColorsSize=2,
    keysize=1, key.title = NA, offsetRow = 0, offsetCol = 0, 
    adjCol = 1, breaks = breaks, symkey = TRUE, labRow = "")
dev.off()
########################################################################################################
# Figure 2B
set.seed(1984)
salmon.counts.pam50 <- cbind(txi.salmon$counts[,rownames(sampleTable)],
  salmon.counts.balancedPrimaries.pam50)
salmon.tpm.pam50 <- cbind(txi.salmon$abundance[,rownames(sampleTable)],
  salmon.tpm.balancedPrimaries.pam50)
number.of.samples.tpmThreshold.pam50 <- apply(salmon.tpm.pam50, MARG = 1, 
  FUN = function(x) length(which(x > 0.5)))
keep.index.pam50 <- as.vector(which(number.of.samples.tpmThreshold.pam50 > 4))

salmon.counts.selected.pam50 <- salmon.counts.pam50[keep.index.pam50,]
y.pam50 <- DGEList(salmon.counts.selected.pam50)
y.pam50 <- calcNormFactors(y.pam50)
cpm.norm.pam50 <- cpm(y.pam50, normalized.lib.sizes = TRUE)
cpm.norm.final.log2.pam50 <- log2(cpm.norm.pam50+1)

sampleTable.pam50 <- rbind(sampleTable[,c(1:6)],sampleTable.balancedPrimaries.pam50)
esr1.pos.df <- subset(sampleTable.pam50, Tumor == "primary" & ER.Status == "Pos")
esr1.neg.df <- subset(sampleTable.pam50, Tumor == "primary" & ER.Status == "Neg")

esr1.pos.patients <- esr1.pos.df$Sample
esr1.neg.patients <- esr1.neg.df$Sample

cpm.norm.final.log2.pam50 <- cpm.norm.final.log2.pam50[pam50.gl,]

ensembl.translation <- getBM(
  filters= "ensembl_gene_id", 
  attributes= c("ensembl_gene_id", "external_gene_name"), 
  values= pam50.gl,
  mart= mart)

row.names(ensembl.translation) <- ensembl.translation$ensembl_gene_id
cpm.norm.final.log2.pam50.merged <- merge(cpm.norm.final.log2.pam50, ensembl.translation,  
  all.x = TRUE, by.x='row.names', by.y= 'ensembl_gene_id')
row.names(cpm.norm.final.log2.pam50.merged) <- cpm.norm.final.log2.pam50.merged$external_gene_name
cpm.norm.final.log2.pam50.merged.final <- cpm.norm.final.log2.pam50.merged[,-c(1,ncol(cpm.norm.final.log2.pam50.merged))]
rownames(cpm.norm.final.log2.pam50.merged.final)[rownames(cpm.norm.final.log2.pam50.merged.final)=="NUF2"] <- "CDCA1"
rownames(cpm.norm.final.log2.pam50.merged.final)[rownames(cpm.norm.final.log2.pam50.merged.final)=="NDC80"] <- "KNTC2"
rownames(cpm.norm.final.log2.pam50.merged.final)[rownames(cpm.norm.final.log2.pam50.merged.final)=="ORC6"] <- "ORC6L"

cpm.pam50 <- t(cpm.norm.final.log2.pam50.merged.final)
# create out df
PAM50.subtype = rep(NA,nrow(cpm.pam50))
pam50.output.df.test <- data.frame(PAM50.subtype.1 = PAM50.subtype, 
  PAM50.subtype.2 = PAM50.subtype,
  PAM50.subtype.3 = PAM50.subtype,
  PAM50.subtype.4 = PAM50.subtype,
  PAM50.subtype.5 = PAM50.subtype,
  PAM50.subtype.6 = PAM50.subtype,
  PAM50.subtype.7 = PAM50.subtype,
  PAM50.subtype.8 = PAM50.subtype,
  PAM50.subtype.9 = PAM50.subtype,
  PAM50.subtype.10 = PAM50.subtype,
  PAM50.subtype.11 = PAM50.subtype,
  PAM50.subtype.12 = PAM50.subtype,
  PAM50.subtype.13 = PAM50.subtype,
  PAM50.subtype.14 = PAM50.subtype,
  PAM50.subtype.15 = PAM50.subtype,
  PAM50.subtype.16 = PAM50.subtype,
  PAM50.subtype.17 = PAM50.subtype,
  PAM50.subtype.18 = PAM50.subtype,
  PAM50.subtype.19 = PAM50.subtype,
  PAM50.subtype.20 = PAM50.subtype,
  row.names = row.names(cpm.pam50))
pam50.probability.df <- data.frame(Basal = rep(NA,43), Her2 = rep(NA,43), 
  LumA = rep(NA,43), LumB = rep(NA,43), Normal = rep(NA,43),  row.names = row.names(cpm.pam50))
pam50.probability.20 <- data.frame(Basal = rep(NA,20), Her2 = rep(NA,20), 
  LumA = rep(NA,20), LumB = rep(NA,20), Normal = rep(NA,20),  row.names = c(1:20))

########## PAM50 CALL
for (i in row.names(cpm.pam50)) {
  for (j in c(1:20)) {
    esr1.pos.patients.test <- esr1.pos.patients[!esr1.pos.patients %in% i]
    esr1.neg.patients.test <- esr1.neg.patients[!esr1.neg.patients %in% i]
      if(identical(esr1.pos.patients.test, esr1.pos.patients)) {
        test.set <- c(sample(esr1.pos.patients.test,7),sample(esr1.neg.patients.test,6))
      }
      if(identical(esr1.neg.patients.test, esr1.neg.patients)) {
        test.set <- c(sample(esr1.pos.patients.test,6),sample(esr1.neg.patients.test,7))
      }
    sample.set <- c(test.set,i)
    pam50.sample.set.final <- cpm.pam50[sample.set,]
    PAM50.subtype <- molecular.subtyping(sbt.model="pam50", 
      data=pam50.sample.set.final, annot=annot.nkis, do.mapping=FALSE)
    pam50.output.df.test[i,j] <- as.character(PAM50.subtype$subtype[i])
    pam50.probability.20[j,] <- PAM50.subtype$subtype.proba[14,]
    print(paste(i, "PAM50 Test", j, "=>",as.character(PAM50.subtype$subtype[i])))
  }
  pam50.probability.df[i,] <- colMeans(pam50.probability.20)
}

pam50.output.df <- as.data.frame(apply(X = pam50.output.df.test, FUN = Mode,MARGIN = 1))
colnames(pam50.output.df) <- "PAM50.subtype"
write.table(pam50.output.df.test, 
  file = "dataOutput/STables/genefu_pam50.bone20foldTest.txt", 
  quote = FALSE, col.names = NA, sep = "\t")
write.table(pam50.output.df, 
  file = "dataOutput/Figure_2/genefu_pam50.boneMet.calls.txt", 
  quote = FALSE, col.names = NA, sep = "\t")
write.table(signif(pam50.probability.df,3), 
  file = "dataOutput/STables/S6.table.txt", 
  quote = FALSE, col.names = NA, sep = "\t")

pam50.subtypes <- data.frame(row.names = rownames(pam50.output.df)[c(1:22)], 
  PAM50.subtype = pam50.output.df[c(1:22),])
pam50.subtypes$Case <- sampleTable$Patient
pam50.subtypes$Tumor <- sampleTable$Tumor
pam50.subtypes$PAM50.plot <- mgsub(c("Basal","Her2", "LumA", "LumB","Normal"), 
  c("1", "2", "3", "4","5"), pam50.subtypes$PAM50.subtype)

pam50.subtypes$Tumor2 <- factor(pam50.subtypes$Tumor, levels = c("primary", "metastasis"))
pam50.subtypes$Case2 <- factor(pam50.subtypes$Case, levels = rev(unique(pam50.subtypes$Case)) )

pdf(file = "dataOutput/Figure_2/pam50-patient-matched-pairs-legend.pdf", width = 5, height = 5 )
plot(c(1,1))
legend("top", legend=c("Basal","Her2","LumA","LumB","Normal"),
  fill=c("#e41a1c","#4daf4a","#377eb8", "#984ea3","goldenrod3"),
  border=FALSE, bty="n", y.intersp = 0.8, x.intersp =0.6, cex=2)
dev.off()

cols2 <- colorRampPalette(c("white", "brown2"))(256)

pam50.probability.df.plot <- pam50.probability.df[rownames(sampleTable),]
pam50.probability.df.plot$sample <- rownames(pam50.probability.df.plot)
pam50.probability.df.melt <- melt(pam50.probability.df.plot, id = "sample")
colnames(pam50.probability.df.melt)[3] <- "Probability"
pam50.probability.df.melt$sample <- factor(pam50.probability.df.melt$sample, 
  levels = rev(rownames(sampleTable)))
p <- ggplot(pam50.probability.df.melt, aes(variable, sample)) + 
  geom_tile(aes(fill = Probability), colour = "white") + 
  scale_fill_gradient(low = "white", high = "brown3") +
  theme(plot.title = element_text(size = 15, face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold"), 
    axis.text.x = element_text(size = 13, face = "bold"),
    axis.title.y = element_text(size=0, face = "bold")) +
  xlab("PAM50 Subtype") + ylab("Sample") + scale_x_discrete(position = "top", expand = c(0,0)) + 
  scale_y_discrete(expand = c(0,0))

pdf(file = "dataOutput/Figure_2/PAM50.probPlot.pdf", height = 10, width = 6)
print(p)
dev.off()

pam50.subtypes$sample <- rownames(pam50.subtypes)
pam50.subtypes$sample <- factor(pam50.subtypes$sample, levels = rev(rownames(sampleTable)))
pdf(file = "dataOutput/Figure_2/pam50.plot.vertical.pdf", width = 1.0, height = 10 )
  qplot(data=pam50.subtypes,x=1,y=sample,fill=factor(PAM50.plot),geom="tile") +
    scale_fill_manual(values=c("1"="#e41a1c", "2"="#4daf4a", "3"="#377eb8", 
      "4" = "#984ea3", "5" = "goldenrod2")) +
    theme_classic() +
    theme(legend.position="none", 
      panel.border=element_blank(),
      plot.title=element_text(lineheight=.8, face="bold"),
      axis.text=element_text(size=5,face="bold"),
      axis.title=element_text(size=4,face="bold"),
      axis.text.x=element_text(size=5,face="bold"),
      axis.title.x=element_text(size=0,face="bold"),
      axis.title.y=element_text(size=0,face="bold")) + 
    scale_y_discrete(expand = c(0,0)) + 
    scale_x_discrete(expand = c(0,0), position = "top") + 
    ggtitle("PAM50")
dev.off()
########################################################################################################
# Figure 2C
cor.results <- rep(NA,11)
names(cor.results) <- unique(sampleTable$Patient)
j <- 1
for (i in seq(1,22,2)) {
  print(i)
  cor.result <- cor(cpm.norm.final.log2[topVarGenes,i],cpm.norm.final.log2[topVarGenes,i+1])
  cor.results[j] <- cor.result
  j <- j + 1
}

bmfs$cor.results <- cor.results
corResult <- cor.test(bmfs$cor.results, bmfs$bmfs)
pdf(file = "dataOutput/Figure_2/bmfs-expressionSimilarity-corr-top10.pdf", width=6, height=6)
plot(cor.results ~ bmfs, data = bmfs, pch = 19, col = color.transparent, 
  main = "Primary/Metastasis Similarity vs. Clinical Time to Metastasis", cex = 2,cex.main = 1, font.lab = 2, 
  ylab = "Expression Correlation Between Patient-Matched Pair", xlab = "Bone Met Free Survival")
abline(fit <- lm(cor.results ~ bmfs, data = bmfs), col='red3', lwd = 3)
legend("topright", bty="n", legend=paste("pearson R = ", 
  signif(as.numeric(corResult$estimate),digits=3), "\n", "p-value < 0.001", sep = ""), text.font = 2)
dev.off()
########################################################################################################
# Differential expression
fdr.rate = 0.10
padjCutoff = 0.10
txi.salmon.final <- txi.salmon
txi.salmon.final$abundance <- txi.salmon$abundance[row.names(cpm.norm.final.log2),sampleTable$Sample]
txi.salmon.final$counts <- txi.salmon$counts[row.names(cpm.norm.final.log2),sampleTable$Sample]
txi.salmon.final$length <- txi.salmon$length[row.names(cpm.norm.final.log2),sampleTable$Sample]

deseq.data <- DESeqDataSetFromTximport(txi = txi.salmon.final,
  colData = sampleTable, design = ~Patient + Tumor)

dds <- DESeq(deseq.data, parallel = TRUE)
genes <- row.names(dds)
results <- results(dds, alpha = fdr.rate, parallel = TRUE, contrast = c("Tumor", "metastasis", "primary" ))
results.padj.ordered <- as.matrix(results[ order(results$padj), ])
geneList <- row.names(results)
ensembl.translation <- getBM(
  filters= "ensembl_gene_id", 
  attributes= c("ensembl_gene_id", "external_gene_name", "external_gene_source", "gene_biotype"), 
  values= geneList,
  mart= mart)
row.names(ensembl.translation) <- ensembl.translation$ensembl_gene_id
results.padj.ordered.merged <- merge(results.padj.ordered, 
  ensembl.translation,  all.x = TRUE, by.x='row.names', by.y= 'ensembl_gene_id')
results.padj.ordered.formatted <- results.padj.ordered.merged[,c(1,8,9,10,2,3,4,5,6,7)]
row.names(results.padj.ordered.formatted) <- results.padj.ordered.formatted$Row.names
results.padj.ordered.final <- results.padj.ordered.formatted[,-1]
save(results.padj.ordered.final, file = "dataOutput/Figure_3/results.padj.ordered.final.Rda")

results.sig <- subset(results.padj.ordered.final, padj < padjCutoff)
results.sig.ordered <- results.sig[ order(results.sig$padj), ]
degenes.10 <- rownames(results.sig)
write(degenes.10, file = "dataOutput/Figure_3/padj-0.10-degenes.txt")

genesUp <- rownames(subset(results.sig.ordered, log2FoldChange > 0))
genesDown <- rownames(subset(results.sig.ordered, log2FoldChange < 0))
DEGs.cpm.norm.log2.up <- cpm.norm.final.log2[genesUp,]
DEGs.cpm.norm.log2.down <- cpm.norm.final.log2[genesDown,]

primary.tumors <- rownames(subset(sampleTable, Tumor == "primary"))
metastatic.tumors <- rownames(subset(sampleTable, Tumor == "metastasis"))
data.hm <- as.matrix(rbind(DEGs.cpm.norm.log2.up[,c(primary.tumors,metastatic.tumors)],
  DEGs.cpm.norm.log2.down[,c(primary.tumors,metastatic.tumors)]))

ensembl.ids.ordered <- ensembl.translation[c(genesUp,genesDown),]
rownames(data.hm) <- ensembl.ids.ordered$external_gene_name

PrimaryMet.colors <- mgsub(c("primary","metastasis"), c("blue2", "firebrick3"), 
  c(rep("primary",11),rep("metastasis",11)))
clab <- cbind(PrimaryMet.colors)
colnames(clab)=c("Tumor")

pdf(file = "dataOutput/Figure_3/DEG-heatmap.pdf", width=6, height=10)
  heatmap.3(data.hm, 
    scale = "row", dendrogram = "none", Colv=FALSE, Rowv=FALSE, ColSideColors=clab,
    col=mycol, cexRow = 0.35, cexCol = 1.0, trace="none", key=FALSE,
    margins=c(10,12), ColSideColorsSize=1,
    keysize=1, key.title = NA, offsetRow = 0, offsetCol = 0, 
    adjCol = 1, breaks = breaks, symkey = TRUE)
dev.off()

# Supplement S7
results.sig.log2FC.ordered <- results.sig[ order(-results.sig$log2FoldChange), ]
genesUp.table <- subset(results.sig.log2FC.ordered, log2FoldChange > 0)
genesDown.table <- subset(results.sig.log2FC.ordered, log2FoldChange < 0)
write.table(genesUp.table, file = "dataOutput/STables/S7.table.genesUp.txt", 
  sep = "\t", quote = FALSE, col.names = NA)
write.table(genesDown.table, file = "dataOutput/STables/S7.table.genesDown.txt", 
  sep = "\t", quote = FALSE, col.names = NA)

# Gene ontology
write(c(genesUp.table$external_gene_name,
  genesDown.table$external_gene_name), 
  file = "dataOutput/STables/S8_table.txt")

########################################################################################################
# GSEA
results.padj.ordered.final.log2FoldChange.ordered <- results.padj.ordered.final[with(results.padj.ordered.final,
  order(-results.padj.ordered.final$log2FoldChange)),]
write.table(results.padj.ordered.final.log2FoldChange.ordered$log2FoldChange,
  file="dataOutput/Figure_4/boneMet.log2FoldChange.rnk",quote=F,sep="\t",
  col.names = FALSE,row.names=results.padj.ordered.final.log2FoldChange.ordered$external_gene_name)
write.table(results.padj.ordered.final.log2FoldChange.ordered$log2FoldChange,
  file="dataOutput/STables/S11.table.txt",quote=F,sep="\t",
  col.names = FALSE,row.names=results.padj.ordered.final.log2FoldChange.ordered$external_gene_name)

# gene plots
# RBBP8 (ENSG00000101773), PIK3C2G (ENSG00000139144), ESR1 (ENSG00000091831),
# FGFR3 (ENSG00000068078), EPHA3 (ENSG00000044524)
pdf(file = "dataOutput/genesOfInterest_pairedPlots.pdf", height = 6, width = 6) 
    pairedPlotList = c("ENSG00000101773", "ENSG00000139144", "ENSG00000091831",
      "ENSG00000068078", "ENSG00000044524")
for (geneOfInterest in pairedPlotList) {
    paired.plot.data <- data.frame(Patient = sampleTable$Patient, 
      Gene = cpm.norm.final.log2[geneOfInterest,], 
      Tumor.Type = sampleTable$Tumor)

    colnames(paired.plot.data)[2] <- geneOfInterest

    paired.plot.data$Tumor.Type <- gsub("primary", "Primary Tumors", 
      paired.plot.data$Tumor.Type)
    paired.plot.data$Tumor.Type <- gsub("metastasis", "Bone Metastases", 
      paired.plot.data$Tumor.Type)
    paired.plot.data$Tumor.Type <- factor(paired.plot.data$Tumor.Type, 
      as.character(paired.plot.data$Tumor.Type))
    paired.plot.data$Tumor.Type <- factor(paired.plot.data$Tumor.Type, 
      c("Primary Tumors", "Bone Metastases"))

    wilcox.result.signedRank <- wilcox.test(paired.plot.data[primary.tumors,2], 
      paired.plot.data[metastatic.tumors,2], paired = TRUE)
    wilcoxSignedRank.table.p.value <- wilcox.result.signedRank$p.value
    wilcoxSignedRank.table.p.value.text <- as.character(signif(wilcoxSignedRank.table.p.value,1))

    print(ggplot(paired.plot.data, aes_string(x="Tumor.Type", 
      y=geneOfInterest, group="Patient", label="Patient")) +
      ggtitle(paste(geneOfInterest, sep = " ")) +
      theme(plot.title = element_text(size = 25, face = "bold"),
        axis.title.x = element_text(size = 0, face = "bold"), 
        axis.text.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size=18, face = "bold"),
        panel.background = element_rect(fill = "gray92")) + 
      geom_point(aes(colour=Tumor.Type), size=7.0, alpha = 0.8, position=position_dodge(width=0.15)) +
      geom_line(size=0.2, alpha=1.0, position=position_dodge(width=0.15)) +
      scale_colour_manual(values=c("blue3", "red3"), guide=FALSE) +
      xlab('Tumor Type') +
      ylab(paste("Log2 Normalized Counts")) + 
      scale_shape_manual(values = c(17, 15,19)) +
      geom_text(x = Inf, y = Inf, label = paste("P=",wilcoxSignedRank.table.p.value.text, sep = " "), 
        hjust = 1, vjust = 1.1, size = 7)
    )
 }   
 dev.off()

# RBBP8 bone-met free survival
survivalSet <- boneMetSurvival.GSE12276
colnames(survivalSet)[3] <- "Met"
expressionSet <- expressionSet.GSE12276.gcrma.medians
geneOfInterest <- "RBBP8"
percentile <- 0.75
maxTime <- 120

brca.gene.expression <- expressionSet[geneOfInterest,]
percentileValue <- quantile(as.numeric(brca.gene.expression), percentile)
brca.gene.expression.high <- which(brca.gene.expression > percentileValue)
brca.gene.expression.low <- which(brca.gene.expression <= percentileValue)
samples.high <- names(expressionSet[brca.gene.expression.high])
samples.low <- names(expressionSet[brca.gene.expression.low])
samples.of.interest.low.high <- c(samples.low, samples.high)
expression.subsets <- data.frame(rep(NA, length(samples.of.interest.low.high)), 
  row.names = samples.of.interest.low.high)
colnames(expression.subsets) <- "expression.subset"
expression.subsets$expression.subset[1:length(samples.low )] <- "low"
expression.subsets$expression.subset[length(samples.low )+1:length(samples.high)] <- "high"
km.data <- transform(merge(survivalSet, expression.subsets, by = 0), 
  row.names=Row.names, Row.names=NULL)
km.data$MFS.months <- as.numeric(km.data$MFS.months)
km.data$Met <- as.numeric(km.data$Met)
km.data$SurvObj <- with(km.data, Surv(pmin(km.data$MFS.months, maxTime), Met == 1))
km.subtype <- survfit(SurvObj ~ expression.subset, data = km.data, conf.type = "log-log")
pdf(file = "dataOutput/Figure_4/RBBP8_bmfs.pdf", width = 8, height = 6)
  print(ggsurvplot(km.subtype,
    pval = TRUE, conf.int = TRUE,
    risk.table = TRUE, 
    risk.table.col = "strata", 
    linetype = 1, 
    surv.median.line = "hv", 
    ggtheme = theme_classic2(), 
    palette = c("firebrick1", "darkblue"),
    risk.table.fontsize = 5.5, risk.table.height = 0.2, 
    surv.plot.height = 0.6, legend = "none", font.x = c(18, "bold", "black"), 
    font.y = c(18, "bold", "black"), font.tickslab = c(15, "plain", "black"), 
    risk.table.title = "# at Risk")
  )
dev.off()
########################################################################################################
# Figure 5
caseIDs <- sampleTable$Patient
salmon.fc <- data.frame(matrix(NA, nrow = nrow(cpm.norm.final.log2), ncol = ncol(cpm.norm.final.log2)/2))
colnames(salmon.fc) <- unique(caseIDs)
rownames(salmon.fc) <- rownames(cpm.norm.final.log2)

for (i in seq(2,ncol(cpm.norm.final.log2),2)) {
  met.sample.id <- colnames(cpm.norm.final.log2)[i]
  primary.sample.id <- colnames(cpm.norm.final.log2)[i-1]
  caseID <- gsub("M", "", met.sample.id)
  salmon.fc[,caseID] <-  cpm.norm.final.log2[,met.sample.id] - cpm.norm.final.log2[,primary.sample.id]
}

ensembl.translation <- getBM(
  filters= "external_gene_name", 
  attributes= c("external_gene_name","ensembl_gene_id"), 
  values= clinicallyActionable.gl,
  mart= mart)

rownames(ensembl.translation) <- ensembl.translation$ensembl_gene_id

genesOfInterest <- intersect(ensembl.translation$ensembl_gene_id, rownames(salmon.fc))
clinicallyActionable.fc <- salmon.fc[genesOfInterest,]
rownames(clinicallyActionable.fc) <- make.unique(ensembl.translation[genesOfInterest,]$external_gene_name)
clinicallyActionable.fc.ordered <- clinicallyActionable.fc[order(rownames(clinicallyActionable.fc)),]
write.table(clinicallyActionable.fc.ordered, file = "dataOutput/STables/STable.13.txt"
  quote = FALSE, col.names = NA, sep = "\t")
########################################################################################################
salmon.fc.melt <- melt(salmon.fc)
colnames(salmon.fc.melt) <- c("Case", "FoldChange")
upRegulated.threshold <- apply(salmon.fc,MARGIN = 2, function(x) quantile(x, probs = 0.95))
downRegulated.threshold <- apply(salmon.fc,MARGIN = 2, function(x) quantile(x, probs = 0.05))
upregulated.matrix <- clinicallyActionable.fc > 0
j <- 1
for (i in upRegulated.threshold) {
  print(i)
  print(j)
  upregulated.matrix[,j] <- clinicallyActionable.fc[,j] > i
  j <- j+1
}
upregulated.matrix.final <- upregulated.matrix[order(-rowSums(upregulated.matrix)),]
class(upregulated.matrix.final) <- "numeric"
downregulated.matrix <- clinicallyActionable.fc < 0
j <- 1
for (i in downRegulated.threshold) {
  print(i)
  print(j)
  downregulated.matrix[,j] <- clinicallyActionable.fc[,j] < i
  j <- j+1
}
downregulated.matrix.final <- downregulated.matrix[order(-rowSums(downregulated.matrix)),]
class(downregulated.matrix.final) <- "numeric"
data.up <- upregulated.matrix.final
data.down <- downregulated.matrix.final
mat.list <- list(Increased = data.up)
mat.list.unified <- unify_mat_list(mat.list)
alter_fun_list = list(
    background = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
    },
    Increased = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "red3", col = NA))
    }
)
col = c("Increased" = "red3")
increased.decreased <- cbind(mat.list.unified$Increased)
recurrent.changes.genes <- names(which(rowSums(increased.decreased) > 3))
mat.list.recurrent <- list(Increased = mat.list.unified$Increased[recurrent.changes.genes,])
pdf("dataOutput/Figure_5/oncoPrint.up.pdf", height = 10, width = 6)
oncoPrint(mat.list.recurrent, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun_list = alter_fun_list, col = col, 
          column_title = "Recurrent Expression Gains in
Clinically Actionable Genes",
          heatmap_legend_param = list(title = "Expression", at = c("Increased"), 
          labels = c("Increased")),row_names_gp = gpar(fontsize = 10,fontface = "bold"), 
          pct_gp = gpar(fontsize =10, fontface = "bold"),
          show_column_names = TRUE)
dev.off()
########### up and down
mat.list <- list(Decreased = data.down)
mat.list.unified <- unify_mat_list(mat.list)
alter_fun_list = list(
    background = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
    },
    Decreased = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "blue3", col = NA))
    }
)
col = c("Decreased" = "blue3")
increased.decreased <- cbind(mat.list.unified$Decreased)
recurrent.changes.genes <- names(which(rowSums(increased.decreased) > 3))
mat.list.recurrent <- list(Decreased = mat.list.unified$Decreased[recurrent.changes.genes,])
pdf("dataOutput/Figure_5/oncoPrint.down.pdf", height = 10, width = 6)
oncoPrint(mat.list.recurrent, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun_list = alter_fun_list, col = col, 
          column_title = "Recurrent Expression Losses in
Clinically Actionable Genes",
          heatmap_legend_param = list(title = "Expression", at = c("Decreased"), 
            labels = c("Decreased")),row_names_gp = gpar(fontsize = 10,fontface = "bold"), 
          pct_gp = gpar(fontsize =10, fontface = "bold"),
          show_column_names = TRUE)
dev.off()
###########
mat.list <- list(Increased = data.up, Decreased = data.down )
mat.list.unified <- unify_mat_list(mat.list)
alter_fun_list = list(
    background = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
    },
    Decreased = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "blue3", col = NA))
    },
    Increased = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "red3", col = NA))
    }
)

col = c("Increased" = "red3", "Decreased" = "blue3")
increased.decreased <- cbind(mat.list.unified$Increased, mat.list.unified$Decreased)
recurrent.changes.genes <- names(which(rowSums(increased.decreased) > 1))
mat.list.recurrent <- list(Increased = mat.list.unified$Increased[recurrent.changes.genes,], 
  Decreased = mat.list.unified$Decreased[recurrent.changes.genes,])
pdf("dataOutput/SFigure_6/oncoPrint.all.pdf", height = 12, width = 5)
oncoPrint(mat.list.recurrent, get_type = function(x) strsplit(x, ";")[[1]],
    alter_fun_list = alter_fun_list, col = col, 
    column_title = "Recurrent Expression Changes in Brain Metastases",
    heatmap_legend_param = list(title = "Expression", at = c("Increased", "Decreased"), 
    labels = c("Increased", "Decreased")),row_names_gp = gpar(fontsize = 5,fontface = "bold"), 
    pct_gp = gpar(fontsize =5, fontface = "bold"), 
    show_column_names = TRUE)
dev.off()

########################################################################################################
caseList <- sampleTable$Patient
plots <- list()
for (singleCase in caseList) {
salmon.fc.melt.case <- subset(salmon.fc.melt, Case == singleCase)
p <- ggplot(salmon.fc.melt.case, aes(FoldChange)) + geom_density(fill = "darkblue",alpha = 0.5) +   
  ggtitle(paste("Case", singleCase, sep = " ")) +
  theme(plot.title = element_text(size = 15, face = "bold"),
    axis.title.x = element_text(size = 10, face = "bold"), 
    axis.text.x = element_text(size = 8, face = "bold"),
    axis.text.y = element_text(size = 8, face = "bold"),
    axis.title.y = element_text(size=10, face = "bold"),
    panel.background = element_rect(fill = "gray92")) +
  geom_vline(xintercept=(downRegulated.threshold[singleCase]), col = "blue3", lwd = 1.0) + 
  geom_vline(xintercept=upRegulated.threshold[singleCase],col = "firebrick3", lwd = 1.0) + 
  xlab("Log2 Fold Change") + ylab("Density") +
  xlim(-4,4)
  plots[[singleCase]] <- p
}
pdf("dataOutput/SFigure_5/caseSpecific-FC-thresholds.pdf", height = 10, width = 10)
  multiplot(plotlist = plots, cols = 3)
dev.off()

########################################################################################################
# end

