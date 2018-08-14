#'@name density
#'@import data.table
#'@export


transcriptplot <- function(DF){ 
  par(mar=c(4,4,4,1))
  plotDensities(log2(DF + 1), legend = "topright", main = "Detected transcripts per sample")
}

