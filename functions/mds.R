#'@name mds
#'@import limma
#'@import data.table
#'@export

mdsplot <- function(DF, DF1, fontsize){
  colnames(DF) <- DF1$shortnames
  plotMDS(as.matrix(DF), dim.plot = c(1,2), cex=fontsize, labels = paste0(DF1$shortnames, ".", DF1$groups), main= "Multi-Dimensional Scaling plot per sample")
}
