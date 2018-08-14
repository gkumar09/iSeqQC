#'@name distribution
#'@import data.table
#'@export

distributionplot <- function(DF){
  par(mar=c(5,10,4,1))
  boxplot(log2(DF + 1), pch=".", horizontal=TRUE, cex.axis=0.9, main= "Distribution of counts per sample", las=1, xlab="log2(counts +1)", col=rainbow(length(unique(DF))))
}
