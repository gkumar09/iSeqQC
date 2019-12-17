#'@name housekeeping
#'@import data.table
#'@export

expression_plot <- function(DF1, DF3){
	gene <- counts_file[grep(paste0("^",input$genes,"$"), counts_file[,1]),]
  if (dim(housekeep)[1] == 2){
    barplot(as.matrix(log2(gene[,-1])), names.arg = as.character(DF1[,2]), las=2, ylim=c(0,50))
    #legend('topright', legend=c("GAPDH", "ACTB"), pch=15, col= c("orange", "blue"))

  }
  else {
    plot(NA, xlim=c(0,1), ylim=c(0,1))
    text(0.5,0.5,"No expression of GAPDH and ACTB found", cex = 1)
}
}