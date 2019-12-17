#'@name phylo
#'@import ape
#'@import data.table
#'@export

phylo_plot <- function(DF3, fontsize){
	plot(as.phylo(hclust(dist(scale(t(DF3)), method="euclidean"))), type="unrooted", main = "Hierarchichal relationship between samples", cex=fontsize)
}
