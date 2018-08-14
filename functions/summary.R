#'@name summary
#'@import psych
#'@import data.table
#'@export


SumMarize<- function(DF){
  sm <- as.data.frame(describe(DF))
  sm=setDT(sm, keep.rownames = T)[]
  sm= sm[,-c(2,6,7,13)]
  colnames(sm) <- c("Samples names", "Number of Genes", "Mean", "SD", "Median", "Min", "Max", "Range", "Skew", "Kurtosis")
  cat("Summary Statistics of the samples\n")
  cat("\n")
  sm
}

distributionplot <- function(DF){
  par(mar=c(5,10,4,1))
  boxplot(log2(DF + 1), pch=".", horizontal=TRUE, cex.axis=0.9, main= "Distribution of counts per sample", las=1, xlab="log2(counts +1)", col=rainbow(length(unique(DF))))
}

transcriptplot <- function(DF){ 
  par(mar=c(4,4,4,1))
  plotDensities(log2(DF + 1), legend = "topright", main = "Detected transcripts per sample")
}

mdsplot <- function(DF, DF1, fontsize){
  colnames(DF) <- DF1$shortnames
  plotMDS(as.matrix(DF), dim.plot = c(1,2), cex=fontsize, labels = paste0(DF1$shortnames, ".", DF1$groups), main= "Multi-Dimensional Scaling plot per sample")
}

pcaplot <- function(DF, DF1){
  respca <- PCA(DF, graph = FALSE)
  fviz_pca_var(respca, col.circle = "white", labelsize= 4, pointsize= 4, label= "var", geom = c("point", "text"), habillage=DF1$groups, title= "Principal Component Variances for all samples") + theme(text = element_text(size = 14))
}

housekeeping_plot <- function(DF1, DF3, DF4, fontsize){
  g1= data.frame(t(DF3[1,]))
  names(g1)<- as.character(DF4[1,1])
  g2= data.frame(t(DF3[2,]))
  names(g2)<- as.character(DF4[2,1])
  g3= data.frame(t(DF3[3,]))
  names(g3)<- as.character(DF4[3,1])
  g4= data.frame(t(DF3[4,]))
  names(g4)<- as.character(DF4[4,1])
  g5= data.frame(t(DF3[5,]))
  names(g5)<- as.character(DF4[5,1])
  g6= data.frame(t(DF3[6,]))
  names(g6)<- as.character(DF4[6,1])
  g7= data.frame(t(DF3[7,]))
  names(g7)<- as.character(DF4[7,1])
  g8= data.frame(t(DF3[8,]))
  names(g8)<- as.character(DF4[8,1])
  g9= data.frame(t(DF3[9,]))
  names(g9)<- as.character(DF4[9,1])
  g10= data.frame(t(DF3[10,]))
  names(g10)<- as.character(DF4[10,1])
  xlim= c(1,dim(DF3)[2])
  plot(NA, xlim= xlim, ylim= c(-5,30), xlab = "Groups", ylab = "CPM", xaxt= "n", main= "CPM variances per sample for the housekeeping genes")
  axis(1, at= 1:dim(DF3)[2], labels = DF1$groups)
  lines(g1, col="#2ca25f")
  lines(g2, col="#8856a7")
  lines(g3, col="#43a2ca")
  lines(g4, col="#e34a33")
  lines(g5, col="#a1d99b")
  lines(g6, col="#000000")
  lines(g7, col="#3f007d")
  lines(g8, col="#f768a1")
  lines(g9, col="#ec7014")
  lines(g10, col="#d0d1e6")
  legend("top", legend = as.character(DF4[,1]), lty = 1, horiz = T,  xpd=T, bty= "n", col= c("#2ca25f", "#8856a7", "#43a2ca", "#e34a33", "#a1d99b", "#000000", "#3f007d", "#f768a1", "#ec7014", "#d0d1e6"), cex = fontsize)
}


phylo_plot <- function(DF3, fontsize){
  plot(as.phylo(hclust(dist(t(DF3)))), type="unrooted", main = "Hierarchichal relationship between samples", cex=fontsize)
}

correlation_plot <- function(DF3){
  cor_cts <- cor(DF3, method="pearson")
  cols<- colorRampPalette(c("red", "white", "blue"))(20)
  corrplot(cor_cts, type="upper", order="hclust", col=cols, method = "color", tl.col="black", tl.srt=45, tl.cex = 0.8, main = "Pearson correlation of all the samples", mar=c(2,2,5,2))
  
}

weight_plot <- function(DF1, DF3){
  par(mar=c(4,8,4,4), las=2)
  x <- voomWithQualityWeights(DF3, plot=F)
  barplot(t(x$targets$sample.weights), col = "grey", horiz = T, names.arg = as.character(DF1$shortnames), main= "Sample specific Weights", cex.names = 0.8)
  abline(v=1, col="red", lty=2)
}