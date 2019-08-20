#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Please provide sample phenotype file and expression matrix file", call.=FALSE)
} else if (length(args)==1) {
  stop("Please provide sample phenotype file and expression matrix file", call.=FALSE)
}

library("shiny")
library("FactoMineR")
library("factoextra")
library("som")
library("psych")
library("data.table")
library("ape")
library("corrplot")
library("limma")

normalize_data <- function(DF){
  data.norm <- som::normalize(DF, byrow=TRUE)
  data.norm <- na.omit(data.norm)
}

filtering_data <- function(DF1, DF2){
  data_new = na.omit(DF2)
  keep <- rowSums(data_new[,-1] >= 5) >= 3
  data_new <- data_new[keep,]
  data_new= data_new[,-1]
  final_data_new=data_new[,grepl(paste0(DF1$samples, collapse = "|"), names(data_new))]
  filt_data_new <- final_data_new[match(as.character(DF1$samples), names(final_data_new))]
}

SumMarize<- function(DF){
  sm <- as.data.frame(describe(DF))
  sm=setDT(sm, keep.rownames = T)[]
  sm= sm[,-c(2,6,7,13)]
  colnames(sm) <- c("Samples names", "Number of Genes", "Mean", "SD", "Median", "Min", "Max", "Range", "Skew", "Kurtosis")
  write.table(sm, file="Basic_Statistics.txt", sep='\t', row.names=F)
}

distributionplot <- function(DF){
	pdf("Counts_distribution.pdf", height = 5, width = 6)
	par(mar=c(5,10,4,1))
	boxplot(log2(DF + 1), pch=".", horizontal=TRUE, cex.axis=0.9, main= "Distribution of counts per sample", las=1, xlab="log2(counts +1)", col="red")
	dev.off()
}

transcriptplot <- function(DF){ 
	geneRpm = t(t(DF[,-1]) * 1e6 / colSums(DF[,-1]))
  	geneRpm = round(geneRpm, 3)
  	pdf("Counts_distribution.pdf", height = 5, width = 6)
  	par(mar=c(4,4,4,1))
  	plotDensities(log2(geneRpm + 0.01), legend = "topright", main = "Mapped Reads density per sample", col=rainbow(length(unique(DF))))
	dev.off()
}

housekeeping_plot <- function(DF1, DF3){
  	housekeepers<- data.frame(gene_symbol= c("GAPDH", "ACTB"))
  	DF3$gene_symbol <- toupper(DF3$gene_symbol)
  	housekeep <- merge(housekeepers, DF3, by="gene_symbol", all=F)
  	pdf("housekeeping_genes.pdf", height = 5, width = 6)
  	barplot(as.matrix(log2(housekeep[,-1])), names.arg = as.character(DF1[,2]), las=2, col=c("orange", "blue"), ylim=c(0,50))
  	legend('topright', legend=c("GADPH", "ACTB"), pch=15, col= c("orange", "blue"))
	dev.off()  
}

phylo_plot <- function(DF3, fontsize){
	pdf("Hierarchical_clustering.pdf", height = 5, width = 6)
  	plot(as.phylo(hclust(dist(t(DF3)))), type="unrooted", main = "Hierarchichal relationship between samples", cex=fontsize)
	dev.off()
}


correlation_plot <- function(DF3){
	cor_cts <- cor(DF3, method="pearson")
  	cols<- colorRampPalette(c("blue", "white", "red"))(20)
  	pdf("Pearson_correlation.pdf", height = 5, width = 6)
  	corrplot(cor_cts, type="upper", order="hclust", col=cols, addCoef.col = "black", method = "color", tl.col="black", tl.srt=45, tl.cex = 0.6, main = "Pearson correlation of all the samples", mar=c(2,2,5,2))
  	dev.off()
}


# basic statistics
manifest_file <- read.delim(args[1], sep="\t", header=TRUE, check.names = F)
manifest_file <- manifest_file[manifest_file$include,]
counts_file <- read.delim(args[2], sep="\t", header=TRUE, check.names = F)
filt_data_new <- filtering_data(manifest_file, counts_file)
colnames(filt_data_new) <- manifest_file$shortnames
data_final <- SumMarize(filt_data_new)


# counts distribution
filt_data_new <- filtering_data(manifest_file, counts_file)
colnames(filt_data_new) <- manifest_file$shortnames
distributionplot(filt_data_new)

# density distribution
filt_data_new <- filtering_data(manifest_file, counts_file)
colnames(filt_data_new) <- manifest_file$shortnames
transcriptplot(filt_data_new)

# housekeeping genes
ncounts_file <- counts_file[,-1]
final_data_new=ncounts_file[,grepl(paste0(manifest_file$samples, collapse = "|"), names(ncounts_file))]
filt_data_new <- final_data_new[match(as.character(manifest_file$samples), names(final_data_new))]
nfilt_data_new <- cbind(counts_file[,1], filt_data_new)
colnames(nfilt_data_new) <- c("gene_symbol", as.character(manifest_file$shortnames))
if (dim(nfilt_data_new)[1]!=2) {
  housekeeping_plot(manifest_file, nfilt_data_new) 
  } else {
  print("Cannot plot housekeeping plots as there is no ACTB and GAPDH housekeeping genes available in your organism")
  } 

# PCA 
filt_data_new <- filtering_data(manifest_file, counts_file) 
data.norm <- normalize_data(filt_data_new)
data.final <- as.matrix(data.norm)
colnames(data.final) <- manifest_file$shortnames
respca <- PCA(data.final, graph = FALSE)
p <- fviz_pca_var(respca, col.circle = rgb(255,255,255, alpha= 0, maxColorValue = 255), labelsize= 4, pointsize= 4, label= "var", geom = c("point", "text"), habillage=manifest_file$groups, title= "Principal Component Variances for all samples") + theme(text = element_text(size = 14))
pdf("PCA_plot.pdf", height = 5, width = 6)
plot(p)
dev.off()

# hierarchiacal clustering
filt_data_new <- filtering_data(manifest_file, counts_file)
data.norm <- normalize_data(filt_data_new)
colnames(data.norm) <- manifest_file$shortnames
phylo_plot(data.norm, 0.9)

# Pearson correlation
filt_data_new <- filtering_data(manifest_file, counts_file)
cor_cts <- cor(filt_data_new, method="pearson")
colnames(cor_cts) <- paste0(manifest_file$shortnames, ".", manifest_file$groups)
row.names(cor_cts) <- paste0(manifest_file$shortnames, ".", manifest_file$groups)
cols<- colorRampPalette(c("blue", "white", "red"))(20)
pdf("Pearson_plot.pdf", height = 5, width = 6)
corrplot(cor_cts, type="upper", order="hclust", col=cols, method = "square",tl.col="black", tl.srt=45, tl.cex = 1, main = "Pearson correlation of all the samples", mar=c(2,2,5,2))
dev.off()

# Spearman correlation
filt_data_new <- filtering_data(manifest_file, counts_file)
cor_cts <- cor(filt_data_new, method="spearman")
colnames(cor_cts) <- paste0(manifest_file$shortnames, ".", manifest_file$groups)
row.names(cor_cts) <- paste0(manifest_file$shortnames, ".", manifest_file$groups)
cols<- colorRampPalette(c("blue", "white", "red"))(20)
pdf("Spearman_plot.pdf", height = 5, width = 6)
corrplot(cor_cts, type="upper", order="hclust", col=cols, tl.col="black", tl.srt=45, tl.cex = 1,method = "square",main = "Spearman correlation of all the samples", mar=c(2,2,5,2))
dev.off()


















