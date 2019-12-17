#'@name filtering
#'@import data.table
#'@export


expfiltering_data <- function(DF1, DF2){
  data_new = na.omit(DF2)
  #keep <- rowSums(data_new[,-1] >= 5) >= 3
  #data_new <- data_new[keep,]
  ndata_new= data_new[,-1]
  final_data_new=ndata_new[,grepl(paste0(DF1$samples, collapse = "|"), names(ndata_new))]
  filt_data_new <- final_data_new[match(as.character(DF1$samples), names(final_data_new))]
  #filt_data_new$gene_symbol <- as.character(data_new[,1])
}
