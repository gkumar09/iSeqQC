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
