#'@name summary
#'@import psych
#'@import data.table
#'@export


SumMarize<- function(DF){
  sm <- as.data.frame(describe(DF, ranges=TRUE, fast = FALSE))
  sm=setDT(sm, keep.rownames = T)[]
  sm= sm[,-c(2,6,7,13)]
  addsum= data.frame(colSums(DF))
  addsum_zero= data.frame(colSums(DF>0))
  newsm= cbind(sm,addsum, addsum_zero)
  colnames(newsm) <- c("Samples", "Detected Genes", "Mean", "SD", "Median", "Min", "Max", "Range", "Skew", "Kurtosis", "Library Size", "Expressed Genes")
  cat("Summary Statistics of the samples\n")
  cat("\n")
  newsm
}
