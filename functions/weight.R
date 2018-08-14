#'@name weight
#'@import limma
#'@import data.table
#'@export

weight_plot <- function(DF1, DF3){
  par(mar=c(4,8,4,4), las=2)
  x <- voomWithQualityWeights(DF3, plot=F)
  barplot(t(x$targets$sample.weights), col = "grey", horiz = T, names.arg = as.character(DF1$shortnames), main= "Sample specific Weights", cex.names = 0.8)
  abline(v=1, col="red", lty=2)
}