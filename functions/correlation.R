#'@name correlation
#'@import corrplot
#'@import data.table
#'@export

correlation_plot <- function(DF3){
  cor_cts <- cor(DF3, method="pearson")
  cols<- colorRampPalette(c("red", "white", "blue"))(20)
  corrplot(cor_cts, type="upper", order="hclust", col=cols, method = "color", tl.col="black", tl.srt=45, tl.cex = 0.8, main = "Pearson correlation of all the samples", mar=c(2,2,5,2))
  
}
