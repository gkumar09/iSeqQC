#'@name housekeeping
#'@import data.table
#'@import shiny
#'@export

housekeeping_plot <- function(DF1, DF3){
	if (input$geneid=="gene_id"){
	  housekeepers<- data.frame(gene_symbol= c("ENSG00000111640", "ENSG00000075624"))
  DF3$gene_symbol <- toupper(DF3$gene_symbol)
  housekeep <- merge(housekeepers, DF3, by="gene_symbol", all=F)
  if (dim(housekeep)[1] == 2){
    barplot(as.matrix(log2(housekeep[,-1])), names.arg = as.character(DF1[,2]), las=2, col=c("orange", "blue"), ylim=c(0,50))
    legend('topright', legend=c("ENSG00000111640", "ENSG00000075624"), pch=15, col= c("orange", "blue"))

  }
  else {
    plot(NA, xlim=c(0,1), ylim=c(0,1))
    text(0.5,0.5,"No expression of GAPDH and ACTB found", cex = 1)
}	
	}
	else if (input$geneid=="gene_symbol"){
	  housekeepers<- data.frame(gene_symbol= c("GAPDH", "ACTB"))
  DF3$gene_symbol <- toupper(DF3$gene_symbol)
  housekeep <- merge(housekeepers, DF3, by="gene_symbol", all=F)
  if (dim(housekeep)[1] == 2){
    barplot(as.matrix(log2(housekeep[,-1])), names.arg = as.character(DF1[,2]), las=2, col=c("orange", "blue"), ylim=c(0,50))
    legend('topright', legend=c("GAPDH", "ACTB"), pch=15, col= c("orange", "blue"))
  }
  else {
    plot(NA, xlim=c(0,1), ylim=c(0,1))
    text(0.5,0.5,"No expression of GAPDH and ACTB found", cex = 1)
}
	
	}

}