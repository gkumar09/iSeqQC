#'@name filtering
#'@import data.table
#'@export


tpmfiltering_data <- function(DF1, DF2){
  data_new = na.omit(DF2)
  data_new= data_new[,-1]
  final_data_new=data_new[,grepl(paste0(DF1$samples, collapse = "|"), names(data_new))]
  filt_data_new <- final_data_new[match(as.character(DF1$samples), names(final_data_new))]
}
