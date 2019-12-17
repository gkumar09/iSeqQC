#'@name normalizing
#'@import som
#'@import data.table
#'@export

normalize_data <- function(DF){
  data.norm <- som::normalize(DF, byrow=TRUE)
  data.norm <- na.omit(data.norm)
}
