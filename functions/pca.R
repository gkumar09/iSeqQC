#'@name pca
#'@import FactoMineR
#'@import factoextra
#'@import data.table
#'@export


pcaplot <- function(DF, DF1){
  respca <- PCA(DF, graph = FALSE)
  fviz_pca_var(respca, col.circle = "white", labelsize= 4, pointsize= 4, label= "var", geom = c("point", "text"), habillage=DF1$groups, title= "Principal Component Variances for all samples") + theme(text = element_text(size = 14))
}
