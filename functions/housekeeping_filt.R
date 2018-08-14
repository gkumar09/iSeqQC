#'@name housekeeping_filt
#'@import data.table
#'@export

housekeeping_filt <- function(DF2, DF3){
  housekeepers<- data.frame(gene_symbol= c("HPRT", "GAPDH", "B2M", "YWHAZ", "TOP1", "GUSB", "RPLP2", "TBP", "PPIA", "ACTB", "UBC"))
  validate(need(names(DF2)[1] == "gene_symbol", "Oops! \n Please re-check the first column name of counts file \n Solution: It should be 'gene_symbol'"))
  housekeep <- merge(housekeepers, DF3, by="gene_symbol")
}