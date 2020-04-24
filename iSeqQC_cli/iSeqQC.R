#!/usr/bin/env Rscript

#########################################
#####
##### usage
##### Rscript --vanilla iSeqQC.R sample_phenotype_file count_matrix type_of_reads type_of_gene_identifier Organism
##### Example: Raw RNA-expression data with gene-symbols as gene identifiers from human
##### Rscript --vanilla iSeqQC.R samplepheno.txt iseqqc_counts.txt R SYMBOL H
##### Example: TPM normalized RNA-expression data with gene-ids as gene identifiers from mouse
##### Rscript --vanilla iSeqQC.R samplepheno.txt iseqqc_counts.txt N ID M
#####
###########################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop("Please provide all necessary arguments i.e. SAMPLE_PHENOTYPE_FILE EXPRESSION_MATRIX FILE_TYPE ORGANISM GENE_IDENTIFIER", call. = FALSE)
}

# library("shiny")
library("FactoMineR")
library("factoextra")
library("som")
library("psych")
library("data.table")
library("ape")
library("corrplot")
library("limma")
# library('EDASeq')  ## Not used

#-------------------------------------------------------------------------------

normalize_data <- function(DF) {
  data.norm <- som::normalize(DF, byrow = TRUE)
  data.norm <- na.omit(data.norm)
}

filtering_data <- function(DF1, DF2) {
  data_new <- na.omit(DF2)
  data_new <- data_new[, -1]
  final_data_new <- data_new[, grepl(paste0(DF1$samples, collapse = "|"), names(data_new))]
  filt_data_new <- final_data_new[match(as.character(DF1$samples), names(final_data_new))]
}

SumMarize <- function(DF) {
  sm <- as.data.frame(describe(DF, ranges = TRUE, fast = FALSE))
  sm <- setDT(sm, keep.rownames = T)[]
  sm <- sm[, -c(2, 6, 7, 13)]
  addsum <- data.frame(colSums(DF))
  addsum_zero <- data.frame(colSums(DF > 0))
  newsm <- cbind(sm, addsum, addsum_zero)
  colnames(newsm) <- c("Samples", "Detected Genes", "Mean", "SD", "Median", "Min", "Max", "Range", "Skew", "Kurtosis", "Library Size", "Expressed Genes")
  write.table(newsm, file = "Summary.Statistics.sample.txt", sep = "\t", row.names = F)
}

distributionplot <- function(DF) {
  pdf("Counts_distribution.pdf", height = 5, width = 6)
  par(mar = c(5, 10, 4, 1))
  boxplot(log2(DF + 1), pch = ".", horizontal = TRUE, cex.axis = 0.9, main = "Distribution of counts per sample", las = 1, xlab = "log2(counts +1)", col = "red")
  dev.off()
}

transcriptplot <- function(DF) {
  geneRpm <- t(t(DF[, -1]) * 1e6 / colSums(DF[, -1]))
  geneRpm <- round(geneRpm, 3)
  pdf("MappedReadDensity_distribution.pdf", height = 5, width = 6)
  par(mar = c(4, 4, 4, 1))
  plotDensities(log2(geneRpm + 0.01), legend = "topright", main = "Mapped Reads density per sample", col = rainbow(length(unique(DF))))
  dev.off()
}

tpmfunction <- function(counts, length) {
  x <- counts / length
  return(t(t(x) * 1e6 / colSums(x, na.rm = T)))
}

phylo_plot <- function(DF3, fontsize) {
  plot(as.phylo(hclust(dist(scale(t(DF3)), method = "euclidean"))), type = "unrooted", main = "Hierarchichal relationship between samples", cex = fontsize)
}

#------------------------------------------------------------------------------

if (args[4] == "ID") {
  gene.col.name <- "gene_id"
  housekeepers <- data.frame(gene_id = c("ENSG00000111640", "ENSG00000075624"))
} else {
  gene.col.name <- "gene_symbol"
  housekeepers <- data.frame(gene_symbol = c("GAPDH", "ACTB"))
}

if (args[5] == "H") {
  nhg <- read.table("datafiles/annotation_hg38withGC.txt", sep = "\t", header = T, check.names = F)
} else if (args[5] == "M") {
  nhg <- read.table("datafiles/annotation_m38withGC.txt", sep = "\t", header = T, check.names = F)
} else {
  nhg <- NULL
}

manifest_file <- read.delim(args[1], sep = "\t", header = TRUE, check.names = F)
manifest_file <- manifest_file[manifest_file$include, ]
counts_file <- read.delim(args[2], sep = "\t", header = TRUE, check.names = F)

#-------------------------------------------------------------------------------

# basic statistics
filt_data_new <- filtering_data(manifest_file, counts_file)
colnames(filt_data_new) <- manifest_file$shortnames
data_final <- SumMarize(filt_data_new)
data_final

# counts distribution
filt_data_new <- filtering_data(manifest_file, counts_file)
colnames(filt_data_new) <- manifest_file$shortnames
distributionplot(filt_data_new)

# mapping reads
filt_data_new <- filtering_data(manifest_file, counts_file)
colnames(filt_data_new) <- manifest_file$shortnames
if (args[4] == "ID") {
  row.names(filt_data_new) <- make.names(gsub("\\..*", "", counts_file[, 1]), unique = T)
  filt_data_new <- setDT(filt_data_new, keep.rownames = T)[]
  names(filt_data_new)[1] <- "gene_id"
} else {
  row.names(filt_data_new) <- make.names(counts_file[, 1], unique = T)
  filt_data_new <- setDT(filt_data_new, keep.rownames = T)[]
  names(filt_data_new)[1] <- "gene_symbol"
}
if (length(nhg)) {
  filt_data_new_withgc <- merge(nhg, filt_data_new, by = gene.col.name)
  final_len <- as.matrix(filt_data_new_withgc$length)
  final_data <- as.data.frame(filt_data_new_withgc[, -c(1, 2, 3, 4)])
  tpm_final_data <- tpmfunction(final_data, final_len)
  pdf("MappedReadDensity_distribution.pdf", height = 5, width = 6)
  par(mar = c(4, 4, 4, 1))
  plotDensities(log2(tpm_final_data + 0.1), legend = "topright", main = "Mapped Reads density per sample", col = rainbow(length(unique(final_data))))
  dev.off()
} else {
  transcriptplot(filt_data_new)
}

# housekeeping
ncounts_file <- counts_file[, -1]
final_data_new <- ncounts_file[, grepl(paste0(manifest_file$samples, collapse = "|"), names(ncounts_file))]
filt_data_new <- final_data_new[match(as.character(manifest_file$samples), names(final_data_new))]
nfilt_data_new <- cbind(counts_file[, 1], filt_data_new)
if (args[4] == "ID") {
  nfilt_data_new[, 1] <- gsub("\\..*", "", nfilt_data_new[, 1])
  colnames(nfilt_data_new) <- c("gene_id", as.character(manifest_file$shortnames))
  housekeepers <- data.frame(gene_id = c("ENSG00000111640", "ENSG00000075624"))
  nfilt_data_new$gene_id <- toupper(nfilt_data_new$gene_id)
  housekeep <- merge(housekeepers, nfilt_data_new, by = "gene_id", all = F)
} else {
  colnames(nfilt_data_new) <- c("gene_symbol", as.character(manifest_file$shortnames))
  housekeepers <- data.frame(gene_symbol = c("GAPDH", "ACTB"))
  nfilt_data_new$gene_symbol <- toupper(nfilt_data_new$gene_symbol)
  housekeep <- merge(housekeepers, nfilt_data_new, by = "gene_symbol", all = F)
}
if (dim(housekeep)[1] == 2) {
  pdf("housekeeping_genes.pdf", height = 5, width = 6)
  barplot(as.matrix(log2(housekeep[, -1])), names.arg = as.character(manifest_file[, 2]), las = 2, col = c("orange", "blue"), ylim = c(0, 50))
  legend("topright", legend = housekeepers[[1]], pch = 15, col = c("orange", "blue"))
  dev.off()
} else {
  pdf("housekeeping_genes.pdf", height = 5, width = 6)
  plot(NA, xlim = c(0, 1), ylim = c(0, 1))
  text(0.5, 0.5, paste0("No expression of ", paste(housekeepers, collapse = " and "), " found"), cex = 1)
  dev.off()
}

# PCA- normalized
filt_data_new <- filtering_data(manifest_file, counts_file)
data.norm <- normalize_data(filt_data_new)
data.final <- as.matrix(data.norm)
colnames(data.final) <- manifest_file$shortnames
respca <- PCA(data.final, graph = FALSE)
p <- fviz_pca_var(respca, col.circle = rgb(255, 255, 255, alpha = 0, maxColorValue = 255), labelsize = 4, pointsize = 4, label = "var", geom = c("point", "text"), habillage = manifest_file$groups, title = "Principal Component Variances for all samples- zscored Normalized") + theme(text = element_text(size = 14))
pdf("PCA_plot_zscored.pdf", height = 5, width = 6)
plot(p)
dev.off()

# PCA- un-normalized
filt_data_new <- filtering_data(manifest_file, counts_file)
data.final <- filt_data_new
colnames(data.final) <- manifest_file$shortnames
respca <- PCA(data.final, graph = FALSE)
p <- fviz_pca_var(respca, col.circle = rgb(255, 255, 255, alpha = 0, maxColorValue = 255), labelsize = 4, pointsize = 4, label = "var", geom = c("point", "text"), habillage = manifest_file$groups, title = "Principal Component Variances for all samples- zscored Normalized") + theme(text = element_text(size = 14))
pdf("PCA_plot_unnormalized.pdf", height = 5, width = 6)
plot(p)
dev.off()

# Multi-factor PCA
filt_data_new <- filtering_data(manifest_file, counts_file)
if (args[3] == "N") {
  data.final <- as.matrix(filt_data_new)
} else {
  data.norm <- normalize_data(filt_data_new)
  data.final <- as.matrix(data.norm)
}
librarysize <- colSums(data.final)
if (ncol(manifest_file) > 4) {
  k <- dim(manifest_file)[2] - 3
  final.tran_data.final <- as.data.frame(cbind(manifest_file[, 5:ncol(manifest_file)], librarysize))
  final.tran_data.final <- apply(final.tran_data.final, 2, as.numeric)
  MFA_respca <- MFA(final.tran_data.final, group = rep(1, k), type = rep("s", k), name.group = c(names(manifest_file)[5:ncol(manifest_file)], "Lib.Size"), graph = FALSE)
  respca <- fviz_mfa_var(MFA_respca, "group", labelsize = 4, pointsize = 4, title = "Multifactor Principal Component Variances") + theme(text = element_text(size = 14))
  pdf("MultifactorPCA_plot.pdf", height = 5, width = 6)
  plot(respca)
  dev.off()
} else {
  pdf(paste0(rpath, "/MultifactorPCA_plot.pdf"), height = 5, width = 6)
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), bty = "n", xaxt = "n", yaxt = "n", xlab = NA, ylab = NA)
  text(0.5, 0.5, "No factor to calculate", cex = 1.2)
  dev.off()
}

# hierarchiacal clustering
filt_data_new <- filtering_data(manifest_file, counts_file)
if (args[3] == "N") {
  data.norm <- filt_data_new
} else {
  data.norm <- normalize_data(filt_data_new)
}
colnames(data.norm) <- manifest_file$shortnames
pdf("Hierarchical_clustering.pdf", height = 5, width = 6)
phylo_plot(data.norm, 0.9)
dev.off()

# Pearson correlation
filt_data_new <- filtering_data(manifest_file, counts_file)
cor_cts <- cor(filt_data_new, method = "pearson")
colnames(cor_cts) <- paste0(manifest_file$shortnames, ".", manifest_file$groups)
row.names(cor_cts) <- paste0(manifest_file$shortnames, ".", manifest_file$groups)
cols <- colorRampPalette(c("blue", "white", "red"))(20)
pdf("Pearson_plot.pdf", height = 5, width = 6)
par(mar = c(2, 2, 5, 2))
corrplot(cor_cts, type = "upper", order = "hclust", col = cols, method = "square", tl.col = "black", tl.srt = 45, tl.cex = 0.7, main = "Pearson correlation of all the samples")
dev.off()

# Spearman correlation
filt_data_new <- filtering_data(manifest_file, counts_file)
cor_cts <- cor(filt_data_new, method = "spearman")
colnames(cor_cts) <- paste0(manifest_file$shortnames, ".", manifest_file$groups)
row.names(cor_cts) <- paste0(manifest_file$shortnames, ".", manifest_file$groups)
cols <- colorRampPalette(c("blue", "white", "red"))(20)
pdf("Spearman_plot.pdf", height = 5, width = 6)
corrplot(cor_cts, type = "upper", order = "hclust", col = cols, tl.col = "black", tl.srt = 45, tl.cex = 0.7, method = "square", main = "Spearman correlation of all the samples", mar = c(2, 2, 5, 2))
dev.off()

# GC-bias
ncts <- na.omit(counts_file)
ndata_new <- ncts[, -1]
final_data_new <- ndata_new[, grepl(paste0(manifest_file$samples, collapse = "|"), names(ndata_new))]
filt_data_new <- final_data_new[match(as.character(manifest_file$samples), names(final_data_new))]
if (args[4] == "ID") {
  row.names(filt_data_new) <- make.names(gsub("\\..*", "", ncts[, 1]), unique = T)
} else {
  row.names(filt_data_new) <- make.names(ncts[, 1], unique = T)
}
log.filt_data_new <- log(filt_data_new)
filt.log.filt_data_new <- log.filt_data_new[!is.infinite(rowSums(log.filt_data_new)), ]
filt.log.filt_data_new <- setDT(filt.log.filt_data_new, keep.rownames = T)[]
names(filt.log.filt_data_new)[1] <- gene.col.name
if (length(nhg)) {
  filt_data_new_withgc <- merge(nhg, filt.log.filt_data_new, by = gene.col.name)
  final_len <- as.matrix(filt_data_new_withgc$length)
  final_gc <- as.matrix(filt_data_new_withgc$gc)
  final_data <- as.matrix(filt_data_new_withgc[, -c(1, 2, 3, 4)])
  cl <- rainbow(dim(final_data)[2])
  pdf("GC_bias.pdf", height = 5, width = 6)
  plot(NA, xlim = c(0.1, 0.9), ylim = c(1, 10), xlab = "GC content", ylab = "Gene counts")
  plotcol <- NA
  for (i in 1:dim(final_data)[2]) {
    points(loess.smooth(final_gc, final_data[, i]), type = "l", col = cl[i], lwd = 1.2)
    plotcol[i] <- cl[i]
  }
  legend("top", legend = manifest_file$shortnames, col = plotcol, lwd = 1.2, cex = 1, bty = "n", xpd = T, ncol = 5)
  dev.off()
} else {
  pdf("GC_bias.pdf", height = 5, width = 6)
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), bty = "n", xaxt = "n", yaxt = "n", xlab = NA, ylab = NA)
  text(0.5, 0.5, "Cannot Plot- the Organism you provided is not compatible", cex = 1.2)
  dev.off()
}
