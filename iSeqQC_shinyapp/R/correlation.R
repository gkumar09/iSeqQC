#'@name correlation
#'@import corrplot
#'@import data.table
#'@export

correlation_plot <- function(DF3){
  cor_cts <- cor(DF3, method="pearson")
  print(cor_cts)
  cols<- colorRampPalette(c("blue", "white", "red"))(20)
  corrplot(cor_cts, type="upper", order="hclust", hclust.method= "ward", col=cols, addCoef.col = "black", method = "color", tl.col="black", tl.srt=45, tl.cex = 0.6, main = "Pearson correlation of all the samples", mar=c(2,2,5,2))
}
