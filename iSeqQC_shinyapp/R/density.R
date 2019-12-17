#'@name density
#'@import limma
#'@export


transcriptplot <- function(DF){ 
  par(mar=c(4,4,4,1))
  genects = apply(DF[,-1],2, function(x) (x/sum(x))*1000000)
  genects = round(genects, 3)
  plotDensities(log2(genects + 0.01), legend = "topright", main = "Mapped Reads density per sample", col=rainbow(length(unique(DF))))
}

tpmfunction <- function(counts,length) {
  x <- counts/length
  return(t(t(x)*1e6/colSums(x, na.rm=T)))
}


