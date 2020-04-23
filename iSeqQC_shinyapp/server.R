library("shiny")
library("FactoMineR")
library("factoextra")
library("som")
library("psych")
library("data.table")
library("ape")
library("corrplot")
library("limma")


options(shiny.maxRequestSize = 1000000*1024^2)


source("R/distribution.R")
source("R/filtering.R")
source("R/normalizing.R")
source("R/phylo.R")
source("R/density.R")
source("R/summary.R")
source("R/housekeeping.R")
source("R/expplot_filtering.R")


shinyServer(function(input, output, session) {


sum_input <- reactive({
    cts_file <- input$countsmatrix
    if (is.null(cts_file)){
      return(invisible())}
    else{
      manifest_file <- read.delim(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
      manifest_file <- manifest_file[manifest_file$include,]
      if (dim(manifest_file)[2]>11){
        return("check manifest file- too many columns")
      }
      else{ 
      counts_file <- read.delim(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
      filt_data_new <- filtering_data(manifest_file, counts_file)
      colnames(filt_data_new) <- manifest_file$shortnames
      data_final <- SumMarize(filt_data_new)
      data_final
      }
    }
}
)

sum2_input <- reactive({
  cts_file <- input$countsmatrix
  if (is.null(cts_file)){
    return(invisible())}
  else{
    manifest_file <- read.delim(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
    manifest_file <- manifest_file[manifest_file$include,]
    counts_file <- read.delim(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
    if (input$ftype=="Raw.Counts"){
      filt_data_new <- filtering_data(manifest_file, counts_file)
    colnames(filt_data_new) <- manifest_file$shortnames
    distributionplot(filt_data_new)}
    else if (input$ftype=="Normalized.Counts"){
      filt_data_new <- filtering_data(manifest_file, counts_file)
      colnames(filt_data_new) <- manifest_file$shortnames
      distributionplot(filt_data_new)
    }
  }
})

sum3_input <- reactive({
  cts_file <- input$countsmatrix
  if (is.null(cts_file)){
    return(invisible())}
  else{
    manifest_file <- read.delim(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
    manifest_file <- manifest_file[manifest_file$include,]
    counts_file <- read.delim(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)

    if (input$ftype=="Raw.Counts"){
      if (input$geneid=="gene_id"){
        ncts = na.omit(counts_file)
        ndata_new= ncts[,-1]
        final_data_new=ndata_new[,grepl(paste0(manifest_file$samples, collapse = "|"), names(ndata_new))]
        filt_data_new <- final_data_new[match(as.character(manifest_file$samples), names(final_data_new))]
        colnames(filt_data_new) <- manifest_file$shortnames
        row.names(filt_data_new) <- make.names(gsub('\\..*','',ncts[,1]), unique=T)
        filt_data_new=setDT(filt_data_new, keep.rownames = T)[]
        names(filt_data_new)[1]<- "gene_id"
          if(input$org=="Human"){
            nhg <- read.table('datafiles/annotation_hg38withGC.txt', sep = '\t', header = T, check.names = F)
          }
          else if(input$org=="Mouse"){
            nhg <- read.table('datafiles/annotation_m38withGC.txt', sep = '\t', header = T, check.names = F)
          }
          else{
            nhg <- data.frame(a=as.Date(character()),
                              b=character(), 
                              c=character(), 
                              stringsAsFactors=FALSE) 
          }
        if (dim(nhg)[1]>1){
          filt_data_new_withgc <- merge(nhg, filt_data_new, by="gene_id")
          final_len <- as.matrix(filt_data_new_withgc$length)
          final_data <- as.data.frame(filt_data_new_withgc[,-c(1,2,3,4)])
          tpm_final_data <- tpmfunction(final_data,final_len)
          plotDensities(log2(tpm_final_data+0.1), legend = "topright", main = "Mapped Reads density per sample", col=rainbow(length(unique(final_data))))
        }
        else {
          transcriptplot(filt_data_new)
        }
        
        }    
      if (input$geneid=="gene_symbol"){
          ncts = na.omit(counts_file)
          ndata_new= ncts[,-1]
          final_data_new=ndata_new[,grepl(paste0(manifest_file$samples, collapse = "|"), names(ndata_new))]
          filt_data_new <- final_data_new[match(as.character(manifest_file$samples), names(final_data_new))]
          colnames(filt_data_new) <- manifest_file$shortnames
          row.names(filt_data_new) <- make.names(ncts[,1], unique=T)
          filt_data_new=setDT(filt_data_new, keep.rownames = T)[]
          names(filt_data_new)[1]<- "gene_symbol"
            if(input$org=="Human"){
              nhg <- read.table('datafiles/annotation_hg38withGC.txt', sep = '\t', header = T, check.names = F)
            }
            else if(input$org=="Mouse"){
              nhg <- read.table('datafiles/annotation_m38withGC.txt', sep = '\t', header = T, check.names = F)
            }
          else{
            nhg <- data.frame(a=as.Date(character()),
                              b=character(), 
                              c=character(), 
                              stringsAsFactors=FALSE) 
          }
          
          if (dim(nhg)[1]>1){
            filt_data_new_withgc <- merge(nhg, filt_data_new, by="gene_symbol")
            final_len <- as.matrix(filt_data_new_withgc$length)
            final_data <- as.data.frame(filt_data_new_withgc[,-c(1,2,3,4)])
            tpm_final_data <- tpmfunction(final_data,final_len)
            plotDensities(log2(tpm_final_data+0.1), legend = "topright", main = "Mapped Reads density per sample", col=rainbow(length(unique(final_data))))
          }
          else {
            transcriptplot(filt_data_new)
          }
       }
    }

    else if (input$ftype=="Normalized.Counts"){
      if (input$geneid=="gene_id"){
        ncts = na.omit(counts_file)
        ndata_new= ncts[,-1]
        final_data_new=ndata_new[,grepl(paste0(manifest_file$samples, collapse = "|"), names(ndata_new))]
        filt_data_new <- final_data_new[match(as.character(manifest_file$samples), names(final_data_new))]
        colnames(filt_data_new) <- manifest_file$shortnames
        row.names(filt_data_new) <- make.names(gsub('\\..*','',ncts[,1]), unique=T)
        filt_data_new=setDT(filt_data_new, keep.rownames = T)[]
        names(filt_data_new)[1]<- "gene_id"
        if(input$org=="Human"){
          nhg <- read.table('datafiles/annotation_hg38withGC.txt', sep = '\t', header = T, check.names = F)
        }
        else if(input$org=="Mouse"){
          nhg <- read.table('datafiles/annotation_m38withGC.txt', sep = '\t', header = T, check.names = F)
        }
        else{
          nhg <- data.frame(a=as.Date(character()),
                            b=character(), 
                            c=character(), 
                            stringsAsFactors=FALSE) 
        }
        
        if (dim(nhg)[1]>1){
          filt_data_new_withgc <- merge(nhg, filt_data_new, by="gene_id")
          final_len <- as.matrix(filt_data_new_withgc$length)
          final_data <- as.data.frame(filt_data_new_withgc[,-c(1,2,3,4)])
          tpm_final_data <- tpmfunction(final_data,final_len)
          plotDensities(log2(tpm_final_data+0.1), legend = "topright", main = "Mapped Reads density per sample", col=rainbow(length(unique(final_data))))
        }
        else {
          transcriptplot(filt_data_new)
        }
        
    }
     
       if (input$geneid=="gene_symbol"){
        ncts = na.omit(counts_file)
        ndata_new= ncts[,-1]
        final_data_new=ndata_new[,grepl(paste0(manifest_file$samples, collapse = "|"), names(ndata_new))]
        filt_data_new <- final_data_new[match(as.character(manifest_file$samples), names(final_data_new))]
        colnames(filt_data_new) <- manifest_file$shortnames
        row.names(filt_data_new) <- make.names(ncts[,1], unique=T)
        filt_data_new=setDT(filt_data_new, keep.rownames = T)[]
        names(filt_data_new)[1]<- "gene_symbol"
        if(input$org=="Human"){
          nhg <- read.table('datafiles/annotation_hg38withGC.txt', sep = '\t', header = T, check.names = F)
        }
        else if(input$org=="Mouse"){
          nhg <- read.table('datafiles/annotation_m38withGC.txt', sep = '\t', header = T, check.names = F)
        }
        else{
          nhg <- data.frame(a=as.Date(character()),
                            b=character(), 
                            c=character(), 
                            stringsAsFactors=FALSE) 
        }
        
        if (dim(nhg)[1]>1){
          filt_data_new_withgc <- merge(nhg, filt_data_new, by="gene_symbol")
          final_len <- as.matrix(filt_data_new_withgc$length)
          final_data <- as.data.frame(filt_data_new_withgc[,-c(1,2,3,4)])
          tpm_final_data <- tpmfunction(final_data,final_len)
          plotDensities(log2(tpm_final_data+0.1), legend = "topright", main = "Mapped Reads density per sample", col=rainbow(length(unique(final_data))))
          
        }
        else {
          transcriptplot(filt_data_new)
        }
        
      }
    }
    }
})

sum4_input <- reactive({
  cts_file <- input$countsmatrix
  if (is.null(cts_file)){
    return(NULL)}
  else{
    manifest_file <- read.table(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
    manifest_file <- manifest_file[manifest_file$include,]
    counts_file <- read.table(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
    ncounts_file <- counts_file[,-1]
    final_data_new=ncounts_file[,grepl(paste0(manifest_file$samples, collapse = "|"), names(ncounts_file))]
    filt_data_new <- final_data_new[match(as.character(manifest_file$samples), names(final_data_new))]
    nfilt_data_new <- cbind(counts_file[,1], filt_data_new)
    if (input$geneid=="gene_id"){
      nfilt_data_new[,1] <- gsub('\\..*','',nfilt_data_new[,1])
      colnames(nfilt_data_new) <- c("gene_id", as.character(manifest_file$shortnames))
      housekeepers<- data.frame(gene_id= c("ENSG00000111640", "ENSG00000075624"))
      nfilt_data_new$gene_id <- toupper(nfilt_data_new$gene_id)
      housekeep <- merge(housekeepers, nfilt_data_new, by="gene_id", all=F)
      if (dim(housekeep)[1] == 2){
        barplot(as.matrix(log2(housekeep[,-1])), names.arg = as.character(manifest_file[,2]), las=2, col=c("orange", "blue"), ylim=c(0,50))
        legend('topright', legend=c("ENSG00000111640", "ENSG00000075624"), pch=15, col= c("orange", "blue"))
        
      }
      else {
        plot(NA, xlim=c(0,1), ylim=c(0,1))
        text(0.5,0.5,"No expression of ENSG00000111640 and ENSG00000075624 found", cex = 1)
      }	
    }
    else if (input$geneid=="gene_symbol"){
      colnames(nfilt_data_new) <- c("gene_symbol", as.character(manifest_file$shortnames))
      housekeepers<- data.frame(gene_symbol= c("GAPDH", "ACTB"))
      nfilt_data_new$gene_symbol <- toupper(nfilt_data_new$gene_symbol)
      housekeep <- merge(housekeepers, nfilt_data_new, by="gene_symbol", all=F)
      if (dim(housekeep)[1] == 2){
        barplot(as.matrix(log2(housekeep[,-1])), names.arg = as.character(manifest_file[,2]), las=2, col=c("orange", "blue"), ylim=c(0,50))
        legend('topright', legend=c("GAPDH", "ACTB"), pch=15, col= c("orange", "blue"))
      }
      else {
        plot(NA, xlim=c(0,1), ylim=c(0,1))
        text(0.5,0.5,"No expression of GAPDH and ACTB found", cex = 1)
      }
    }
  }
})

pca2_input <- reactive({
  cts_file <- input$countsmatrix
  if (is.null(cts_file)){
    return(NULL)}
  else{
    manifest_file <- read.delim(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
    manifest_file <- manifest_file[manifest_file$include,]
    counts_file <- read.delim(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
    if (input$ftype=="Raw.Counts"){
    filt_data_new <- filtering_data(manifest_file, counts_file) 
    data.norm <- normalize_data(filt_data_new)
    data.final <- as.matrix(data.norm)
    colnames(data.final) <- manifest_file$shortnames
    respca <- PCA(data.final, graph = FALSE)
    p <- fviz_pca_var(respca, col.circle = rgb(255,255,255, alpha= 0, maxColorValue = 255), labelsize= 4, pointsize= 4, label= "var", geom = c("point", "text"), habillage=manifest_file$groups, title= "Principal Component Variances for all samples- zscored Normalized") + theme(text = element_text(size = 14))
    plot(p)
    }
    else if (input$ftype=="Normalized.Counts"){
      filt_data_new <- filtering_data(manifest_file, counts_file)
      data.norm <- normalize_data(filt_data_new)
      data.final <- as.matrix(data.norm)
      data.final <- filt_data_new
      colnames(data.final) <- manifest_file$shortnames
      respca <- PCA(data.final, graph = FALSE)
      p <- fviz_pca_var(respca, col.circle = rgb(255,255,255, alpha= 0, maxColorValue = 255), labelsize= 4, pointsize= 4, label= "var", geom = c("point", "text"), habillage=manifest_file$groups, title= "Principal Component Variances for all samples- zscored Normalized") + theme(text = element_text(size = 14))
      plot(p)
    }
      
      }
})

pca3_input <- reactive({
  cts_file <- input$countsmatrix
  if (is.null(cts_file)){
    return(NULL)}
  else{
    manifest_file <- read.delim(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
    manifest_file <- manifest_file[manifest_file$include,]
    counts_file <- read.delim(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
    if (input$ftype=="Raw.Counts"){
      filt_data_new <- filtering_data(manifest_file, counts_file) 
      data.final <- filt_data_new
      colnames(data.final) <- manifest_file$shortnames
      respca <- PCA(data.final, graph = FALSE)
      p <- fviz_pca_var(respca, col.circle = rgb(255,255,255, alpha= 0, maxColorValue = 255), labelsize= 4, pointsize= 4, label= "var", geom = c("point", "text"), habillage=manifest_file$groups, title= "Principal Component Variances for all samples- Un-Normalized") + theme(text = element_text(size = 14))
      plot(p)
    }
    else if (input$ftype=="Normalized.Counts"){
      filt_data_new <- filtering_data(manifest_file, counts_file) 
      data.final <- filt_data_new
      colnames(data.final) <- manifest_file$shortnames
      respca <- PCA(data.final, graph = FALSE)
      p <- fviz_pca_var(respca, col.circle = rgb(255,255,255, alpha= 0, maxColorValue = 255), labelsize= 4, pointsize= 4, label= "var", geom = c("point", "text"), habillage=manifest_file$groups, title= "Principal Component Variances for all samples- Un-Normalized") + theme(text = element_text(size = 14))
      plot(p)
    }
    }
})

pca4_input <- reactive({
  cts_file <- input$countsmatrix
  if (is.null(cts_file)){
    return(NULL)}
  else{
    manifest_file <- read.delim(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
    manifest_file <- manifest_file[manifest_file$include,]
    counts_file <- read.delim(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
    filt_data_new <- filtering_data(manifest_file, counts_file) 
	data.norm <- normalize_data(filt_data_new)
	data.final <- as.matrix(data.norm)
	librarysize <- colSums(data.final)
      if(dim(manifest_file)[2]==5){
        final.tran_data.final <- as.data.frame(cbind(as.numeric(manifest_file[,5]),librarysize))
        MFA_respca <- MFA(final.tran_data.final, group = c(1,1), type = c('s','s'), name.group = c(names(manifest_file[5]),"Lib.Size"), graph = F) 
        respca <- fviz_mfa_var(MFA_respca, "group", labelsize= 3, pointsize= 4, title= "Multifactor Principal Component Variances")+ theme(text = element_text(size = 14))
        plot(respca)
        
      }
      else if(dim(manifest_file)[2]==6){
        final.tran_data.final <- as.data.frame(cbind(as.numeric(manifest_file[,5]),as.numeric(manifest_file[,6]),librarysize))
        MFA_respca <- MFA(final.tran_data.final, group = c(1,1,1), type = c('s','s','s'), name.group = c(names(manifest_file[5]),names(manifest_file[6]),"Lib.Size"), graph = F) 
        respca <- fviz_mfa_var(MFA_respca, "group", labelsize= 4, pointsize= 4, title= "Multifactor Principal Component Variances")+ theme(text = element_text(size = 14))
        plot(respca)
      }
      else if(dim(manifest_file)[2]==7){
        final.tran_data.final <- as.data.frame(cbind(as.numeric(manifest_file[,5]),as.numeric(manifest_file[,6]),as.numeric(manifest_file[,7]),librarysize))
        MFA_respca <- MFA(final.tran_data.final, group = c(1,1,1,1), type = c('s','s','s','s'), name.group = c(names(manifest_file[5]),names(manifest_file[6]),names(manifest_file[7]),"Lib.Size"), graph = F) 
        respca <- fviz_mfa_var(MFA_respca, "group", labelsize= 4, pointsize= 4, title= "Multifactor Principal Component Variances")+ theme(text = element_text(size = 14))
        plot(respca)
      }
      else if(dim(manifest_file)[2]==8){
        final.tran_data.final <- as.data.frame(cbind(as.numeric(manifest_file[,5]),as.numeric(manifest_file[,6]),as.numeric(manifest_file[,7]),as.numeric(manifest_file[,8]),librarysize))
        MFA_respca <- MFA(final.tran_data.final, group = c(1,1,1,1,1), type = c('s','s','s','s','s'), name.group = c(names(manifest_file[5]),names(manifest_file[6]),names(manifest_file[7]),names(manifest_file[8]),"Lib.Size"), graph = F) 
        respca <- fviz_mfa_var(MFA_respca, "group", labelsize= 4, pointsize= 4, title= "Multifactor Principal Component Variances")+ theme(text = element_text(size = 14))
        plot(respca)
      }
      else if(dim(manifest_file)[2]==9){
        final.tran_data.final <- as.data.frame(cbind(as.numeric(manifest_file[,5]),as.numeric(manifest_file[,6]),as.numeric(manifest_file[,7]),as.numeric(manifest_file[,8]),as.numeric(manifest_file[,9]),librarysize))
        MFA_respca <- MFA(final.tran_data.final, group = c(1,1,1,1,1,1), type = c('s','s','s','s','s','s'), name.group = c(names(manifest_file[5]),names(manifest_file[6]),names(manifest_file[7]),names(manifest_file[8]),names(manifest_file[9]),"Lib.Size"), graph = F) 
        respca <- fviz_mfa_var(MFA_respca, "group", labelsize= 4, pointsize= 4, title= "Multifactor Principal Component Variances")+ theme(text = element_text(size = 14))
        plot(respca)
      }
      else if(dim(manifest_file)[2]==10){
        final.tran_data.final <- as.data.frame(cbind(as.numeric(manifest_file[,5]),as.numeric(manifest_file[,6]),as.numeric(manifest_file[,7]),as.numeric(manifest_file[,8]),as.numeric(manifest_file[,9]),as.numeric(manifest_file[,10]),librarysize))
        MFA_respca <- MFA(final.tran_data.final, group = c(1,1,1,1,1,1,1), type = c('s','s','s','s','s','s','s'), name.group = c(names(manifest_file[5]),names(manifest_file[6]),names(manifest_file[7]),names(manifest_file[8]),names(manifest_file[9]),names(manifest_file[10]),"Lib.Size"), graph = F) 
        respca <- fviz_mfa_var(MFA_respca, "group", labelsize= 4, pointsize= 4, title= "Multifactor Principal Component Variances")+ theme(text = element_text(size = 14))
        plot(respca)
      }
      else if(dim(manifest_file)[2]==11){
        final.tran_data.final <- as.data.frame(cbind(as.numeric(manifest_file[,5]),as.numeric(manifest_file[,6]),as.numeric(manifest_file[,7]),as.numeric(manifest_file[,8]),as.numeric(manifest_file[,9]),as.numeric(manifest_file[,10]),as.numeric(manifest_file[,11]),librarysize))
        MFA_respca <- MFA(final.tran_data.final, group = c(1,1,1,1,1,1,1,1), type = c('s','s','s','s','s','s','s','s'), name.group = c(names(manifest_file[5]),names(manifest_file[6]),names(manifest_file[7]),names(manifest_file[8]),names(manifest_file[9]),names(manifest_file[10]),names(manifest_file[11]),"Lib.Size"), graph = F) 
        respca <- fviz_mfa_var(MFA_respca, "group", labelsize= 4, pointsize= 4, title= "Multifactor Principal Component Variances")+ theme(text = element_text(size = 14))
        plot(respca)
      }
      else{
        plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
        text(0.5,0.5,"No factor to calculate", cex = 1.2)
       }
    }
})


cor1_input <- reactive({
  cts_file <- input$countsmatrix
  if (is.null(cts_file)){
    return(NULL)}
  else{
    manifest_file <- read.delim(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
    manifest_file <- manifest_file[manifest_file$include,]
    counts_file <- read.delim(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
    if (input$ftype=="Raw.Counts"){
    filt_data_new <- filtering_data(manifest_file, counts_file)
    filt_data_new <- normalize_data(filt_data_new)
    colnames(filt_data_new) <- manifest_file$shortnames
    phylo_plot(filt_data_new, 0.9)
    }
    else if (input$ftype=="Normalized.Counts"){
      filt_data_new <- filtering_data(manifest_file, counts_file)
      colnames(filt_data_new) <- manifest_file$shortnames
      phylo_plot(filt_data_new, 0.9)
    }
    
    }
})

cor2_input <- reactive({
  cts_file <- input$countsmatrix
  if (is.null(cts_file)){
    return(NULL)}
  else{
    manifest_file <- read.delim(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
    manifest_file <- manifest_file[manifest_file$include,]
    counts_file <- read.delim(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
    if (input$ftype=="Raw.Counts"){
    filt_data_new <- filtering_data(manifest_file, counts_file)
    cor_cts <- cor(filt_data_new, method="pearson")
    colnames(cor_cts) <- paste0(manifest_file$shortnames, ".", manifest_file$groups)
    row.names(cor_cts) <- paste0(manifest_file$shortnames, ".", manifest_file$groups)
    cols<- colorRampPalette(c("blue", "white", "red"))(20)
    corrplot(cor_cts, type="upper", order="hclust", col=cols, method = "square",tl.col="black", tl.srt=45, tl.cex = 1, main = "Pearson correlation of all the samples", mar=c(2,2,5,2))
    }
    else if (input$ftype=="Normalized.Counts"){
      filt_data_new <- filtering_data(manifest_file, counts_file)
      cor_cts <- cor(filt_data_new, method="pearson")
      colnames(cor_cts) <- paste0(manifest_file$shortnames, ".", manifest_file$groups)
      row.names(cor_cts) <- paste0(manifest_file$shortnames, ".", manifest_file$groups)
      cols<- colorRampPalette(c("blue", "white", "red"))(20)
      corrplot(cor_cts, type="upper", order="hclust", col=cols, method = "square",tl.col="black", tl.srt=45, tl.cex = 1, main = "Pearson correlation of all the samples", mar=c(2,2,5,2))
    }
    
    }
})

cor3_input <- reactive({
  cts_file <- input$countsmatrix
  if (is.null(cts_file)){
    return(NULL)}
  else{
    manifest_file <- read.delim(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
    manifest_file <- manifest_file[manifest_file$include,]
    counts_file <- read.delim(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
    if (input$ftype=="Raw.Counts"){
    filt_data_new <- filtering_data(manifest_file, counts_file)
    cor_cts <- cor(filt_data_new, method="spearman")
    colnames(cor_cts) <- paste0(manifest_file$shortnames, ".", manifest_file$groups)
    row.names(cor_cts) <- paste0(manifest_file$shortnames, ".", manifest_file$groups)
    cols<- colorRampPalette(c("blue", "white", "red"))(20)
    corrplot(cor_cts, type="upper", order="hclust", col=cols, tl.col="black", tl.srt=45, tl.cex = 1,method = "square",main = "Spearman correlation of all the samples", mar=c(2,2,5,2))
    }
    else if (input$ftype=="Normalized.Counts"){
      filt_data_new <- filtering_data(manifest_file, counts_file)
      cor_cts <- cor(filt_data_new, method="spearman")
      colnames(cor_cts) <- paste0(manifest_file$shortnames, ".", manifest_file$groups)
      row.names(cor_cts) <- paste0(manifest_file$shortnames, ".", manifest_file$groups)
      cols<- colorRampPalette(c("blue", "white", "red"))(20)
      corrplot(cor_cts, type="upper", order="hclust", col=cols, tl.col="black", tl.srt=45, tl.cex = 1,method = "square",main = "Spearman correlation of all the samples", mar=c(2,2,5,2))
    }
    
    }
})

exp_input <- reactive({
  cts_file <- input$countsmatrix
  if (is.null(cts_file)){
    return(NULL)}
  else{
    manifest_file <- read.table(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
    manifest_file <- manifest_file[manifest_file$include,]
    counts_file <- read.table(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
    data_new = na.omit(counts_file)
    ndata_new= data_new[,-1]
    final_data_new=ndata_new[,grepl(paste0(manifest_file$samples, collapse = "|"), names(ndata_new))]
    filt_data_new <- final_data_new[match(as.character(manifest_file$samples), names(final_data_new))]
    data_new[,1] <- gsub('\\..*','',data_new[,1])
    row.names(filt_data_new) <- make.names(as.character(data_new[,1]), unique=T)
    colnames(filt_data_new) <- manifest_file$shortnames
    gene_query <- input$genes
    gene <- filt_data_new[grep(paste0("^",toupper(input$genes),"$"), toupper(row.names(filt_data_new))),]
    if (dim(gene)[1] == 1){
      barplot(as.matrix(log2(gene)), names.arg = as.character(manifest_file$shortnames), main= as.character(input$genes),las=2, ylab="Log2 expression", col="#fb6a4a", ylim=c(0,max(log2(gene[,-1]+5))))
    }
    else if (dim(gene)[1] == 0){
      plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
      text(0.5,0.5,"Please choose gene IDs", cex = 1.2)
    }
    else {
      plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
      text(0.5,0.5,"Cannot Plot- check if you are using right gene identifier ", cex = 1.2)
    }
  }
})

bias_input <- reactive({
  cts_file <- input$countsmatrix
  if (is.null(cts_file)){
    return(NULL)}
  else{
    manifest_file <- read.table(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
    manifest_file <- manifest_file[manifest_file$include,]
    counts_file <- read.table(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
    ncts = na.omit(counts_file)
    ndata_new= ncts[,-1]
    final_data_new=ndata_new[,grepl(paste0(manifest_file$samples, collapse = "|"), names(ndata_new))]
    filt_data_new <- final_data_new[match(as.character(manifest_file$samples), names(final_data_new))]
    
    if(input$geneid=="gene_id"){
      row.names(filt_data_new) <- make.names(gsub('\\..*','',ncts[,1]), unique=T)
      log.filt_data_new <- log(filt_data_new)
      filt.log.filt_data_new <- log.filt_data_new[!is.infinite(rowSums(log.filt_data_new)),]
      filt.log.filt_data_new =setDT(filt.log.filt_data_new, keep.rownames = T)[]
      names(filt.log.filt_data_new)[1]<- "gene_id"
      
      if(input$org=="Human"){
        nhg <- read.table('datafiles/annotation_hg38withGC.txt', sep = '\t', header = T, check.names = F)
      }
      else if(input$org=="Mouse"){
        nhg <- read.table('datafiles/annotation_m38withGC.txt', sep = '\t', header = T, check.names = F)
      }
      
      else{
        nhg <- data.frame(a=as.Date(character()),
                          b=character(), 
                          c=character(), 
                          stringsAsFactors=FALSE) 
      }
      
      if (dim(nhg)[1]>1){
        filt_data_new_withgc <- merge(nhg, filt.log.filt_data_new, by="gene_id")
        final_len <- as.matrix(filt_data_new_withgc$length)
        final_gc <- as.matrix(filt_data_new_withgc$gc)
        final_data <- as.matrix(filt_data_new_withgc[,-c(1,2,3,4)])
        cl <- rainbow(dim(final_data)[2])
        plot(NA, xlim=c(0.1,0.9),ylim=c(1,10), xlab="GC content", ylab="Gene counts")
        plotcol=NA
        for (i in 1:dim(final_data)[2]){
          points(loess.smooth(final_gc,final_data[,i]), type='l', col=cl[i],lwd=1.2)
          plotcol[i] <- cl[i]
        }
        legend("top", legend = manifest_file$shortnames, col = plotcol, lwd = 1.2, cex = 1, bty='n', xpd=T, ncol=5)
      }
      else {
        plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
        text(0.5,0.5,"Cannot Plot- the Organism you provided is not compatible", cex = 1.2)
      }
      }
      

      else if(input$geneid=="gene_symbol"){
        row.names(filt_data_new) <- make.names(ncts[,1], unique=T)
        log.filt_data_new <- log(filt_data_new)
        filt.log.filt_data_new <- log.filt_data_new[!is.infinite(rowSums(log.filt_data_new)),]
        filt.log.filt_data_new =setDT(filt.log.filt_data_new, keep.rownames = T)[]
        names(filt.log.filt_data_new)[1]<- "gene_symbol"
        
        if(input$org=="Human"){
          nhg <- read.table('datafiles/annotation_hg38withGC.txt', sep = '\t', header = T, check.names = F)
        }
        else if(input$org=="Mouse"){
          nhg <- read.table('datafiles/annotation_m38withGC.txt', sep = '\t', header = T, check.names = F)
        }
        else{
          nhg <- data.frame(a=as.Date(character()),
                            b=character(), 
                            c=character(), 
                            stringsAsFactors=FALSE) 
        }
        
        if (dim(nhg)[1]>1){
          
          filt_data_new_withgc <- merge(nhg, filt.log.filt_data_new, by="gene_symbol")
          final_len <- as.matrix(filt_data_new_withgc$length)
          final_gc <- as.matrix(filt_data_new_withgc$gc)
          final_data <- as.matrix(filt_data_new_withgc[,-c(1,2,3,4)])
          cl <- rainbow(dim(final_data)[2])
          plot(NA, xlim=c(0.1,0.9),ylim=c(1,10), xlab="GC content", ylab="Gene counts")
          plotcol=NA
          for (i in 1:dim(final_data)[2]){
            points(loess.smooth(final_gc,final_data[,i]), type='l', col=cl[i],lwd=1.2)
            plotcol[i] <- cl[i]
          }
          legend("top", legend = manifest_file$shortnames, col = plotcol, lwd = 1.2, cex = 1, bty='n', xpd=T, ncol=5)
        }
        
        else {
          plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
          text(0.5,0.5,"Cannot Plot- the Organism you provided is not compatible", cex = 1.2)
        }
   
    }
  }
})


# downloading plots
select_plot1 <- function(){
    manifest_file <- read.delim(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
    manifest_file <- manifest_file[manifest_file$include,]
    counts_file <- read.delim(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
    if (input$ftype=="Raw.Counts"){
      filt_data_new <- filtering_data(manifest_file, counts_file)
    colnames(filt_data_new) <- manifest_file$shortnames
    distributionplot(filt_data_new)}
    else if (input$ftype=="Normalized.Counts"){
      filt_data_new <- filtering_data(manifest_file, counts_file)
      colnames(filt_data_new) <- manifest_file$shortnames
      distributionplot(filt_data_new)
  }
}


select_plot2 <- function(){
    manifest_file <- read.delim(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
    manifest_file <- manifest_file[manifest_file$include,]
    counts_file <- read.delim(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)

    if (input$ftype=="Raw.Counts"){
      if (input$geneid=="gene_id"){
        ncts = na.omit(counts_file)
        ndata_new= ncts[,-1]
        final_data_new=ndata_new[,grepl(paste0(manifest_file$samples, collapse = "|"), names(ndata_new))]
        filt_data_new <- final_data_new[match(as.character(manifest_file$samples), names(final_data_new))]
        colnames(filt_data_new) <- manifest_file$shortnames
        row.names(filt_data_new) <- make.names(gsub('\\..*','',ncts[,1]), unique=T)
        filt_data_new=setDT(filt_data_new, keep.rownames = T)[]
        names(filt_data_new)[1]<- "gene_id"
          if(input$org=="Human"){
            nhg <- read.table('datafiles/annotation_hg38withGC.txt', sep = '\t', header = T, check.names = F)
          }
          else if(input$org=="Mouse"){
            nhg <- read.table('datafiles/annotation_m38withGC.txt', sep = '\t', header = T, check.names = F)
          }
          else{
            nhg <- data.frame(a=as.Date(character()),
                              b=character(), 
                              c=character(), 
                              stringsAsFactors=FALSE) 
          }
        if (dim(nhg)[1]>1){
          filt_data_new_withgc <- merge(nhg, filt_data_new, by="gene_id")
          final_len <- as.matrix(filt_data_new_withgc$length)
          final_data <- as.data.frame(filt_data_new_withgc[,-c(1,2,3,4)])
          tpm_final_data <- tpmfunction(final_data,final_len)
          plotDensities(log2(tpm_final_data+0.1), legend = "topright", main = "Mapped Reads density per sample", col=rainbow(length(unique(final_data))))
        }
        else {
          transcriptplot(filt_data_new)
        }
        
        }    
      if (input$geneid=="gene_symbol"){
          ncts = na.omit(counts_file)
          ndata_new= ncts[,-1]
          final_data_new=ndata_new[,grepl(paste0(manifest_file$samples, collapse = "|"), names(ndata_new))]
          filt_data_new <- final_data_new[match(as.character(manifest_file$samples), names(final_data_new))]
          colnames(filt_data_new) <- manifest_file$shortnames
          row.names(filt_data_new) <- make.names(ncts[,1], unique=T)
          filt_data_new=setDT(filt_data_new, keep.rownames = T)[]
          names(filt_data_new)[1]<- "gene_symbol"
            if(input$org=="Human"){
              nhg <- read.table('datafiles/annotation_hg38withGC.txt', sep = '\t', header = T, check.names = F)
            }
            else if(input$org=="Mouse"){
              nhg <- read.table('datafiles/annotation_m38withGC.txt', sep = '\t', header = T, check.names = F)
            }
          else{
            nhg <- data.frame(a=as.Date(character()),
                              b=character(), 
                              c=character(), 
                              stringsAsFactors=FALSE) 
          }
          
          if (dim(nhg)[1]>1){
            filt_data_new_withgc <- merge(nhg, filt_data_new, by="gene_symbol")
            final_len <- as.matrix(filt_data_new_withgc$length)
            final_data <- as.data.frame(filt_data_new_withgc[,-c(1,2,3,4)])
            tpm_final_data <- tpmfunction(final_data,final_len)
            plotDensities(log2(tpm_final_data+0.1), legend = "topright", main = "Mapped Reads density per sample", col=rainbow(length(unique(final_data))))
          }
          else {
            transcriptplot(filt_data_new)
          }
       }
    }

    else if (input$ftype=="Normalized.Counts"){
      if (input$geneid=="gene_id"){
        ncts = na.omit(counts_file)
        ndata_new= ncts[,-1]
        final_data_new=ndata_new[,grepl(paste0(manifest_file$samples, collapse = "|"), names(ndata_new))]
        filt_data_new <- final_data_new[match(as.character(manifest_file$samples), names(final_data_new))]
        colnames(filt_data_new) <- manifest_file$shortnames
        row.names(filt_data_new) <- make.names(gsub('\\..*','',ncts[,1]), unique=T)
        filt_data_new=setDT(filt_data_new, keep.rownames = T)[]
        names(filt_data_new)[1]<- "gene_id"
        if(input$org=="Human"){
          nhg <- read.table('datafiles/annotation_hg38withGC.txt', sep = '\t', header = T, check.names = F)
        }
        else if(input$org=="Mouse"){
          nhg <- read.table('datafiles/annotation_m38withGC.txt', sep = '\t', header = T, check.names = F)
        }
        else{
          nhg <- data.frame(a=as.Date(character()),
                            b=character(), 
                            c=character(), 
                            stringsAsFactors=FALSE) 
        }
        
        if (dim(nhg)[1]>1){
          filt_data_new_withgc <- merge(nhg, filt_data_new, by="gene_id")
          final_len <- as.matrix(filt_data_new_withgc$length)
          final_data <- as.data.frame(filt_data_new_withgc[,-c(1,2,3,4)])
          tpm_final_data <- tpmfunction(final_data,final_len)
          plotDensities(log2(tpm_final_data+0.1), legend = "topright", main = "Mapped Reads density per sample", col=rainbow(length(unique(final_data))))
        }
        else {
          transcriptplot(filt_data_new)
        }
        
    }
     
       if (input$geneid=="gene_symbol"){
        ncts = na.omit(counts_file)
        ndata_new= ncts[,-1]
        final_data_new=ndata_new[,grepl(paste0(manifest_file$samples, collapse = "|"), names(ndata_new))]
        filt_data_new <- final_data_new[match(as.character(manifest_file$samples), names(final_data_new))]
        colnames(filt_data_new) <- manifest_file$shortnames
        row.names(filt_data_new) <- make.names(ncts[,1], unique=T)
        filt_data_new=setDT(filt_data_new, keep.rownames = T)[]
        names(filt_data_new)[1]<- "gene_symbol"
        if(input$org=="Human"){
          nhg <- read.table('datafiles/annotation_hg38withGC.txt', sep = '\t', header = T, check.names = F)
        }
        else if(input$org=="Mouse"){
          nhg <- read.table('datafiles/annotation_m38withGC.txt', sep = '\t', header = T, check.names = F)
        }
        else{
          nhg <- data.frame(a=as.Date(character()),
                            b=character(), 
                            c=character(), 
                            stringsAsFactors=FALSE) 
        }
        
        if (dim(nhg)[1]>1){
          filt_data_new_withgc <- merge(nhg, filt_data_new, by="gene_symbol")
          final_len <- as.matrix(filt_data_new_withgc$length)
          final_data <- as.data.frame(filt_data_new_withgc[,-c(1,2,3,4)])
          tpm_final_data <- tpmfunction(final_data,final_len)
          plotDensities(log2(tpm_final_data+0.1), legend = "topright", main = "Mapped Reads density per sample", col=rainbow(length(unique(final_data))))
          
        }
        else {
          transcriptplot(filt_data_new)
        }

    }
  }
}


select_plot3 <- function(){
    manifest_file <- read.table(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
    manifest_file <- manifest_file[manifest_file$include,]
    counts_file <- read.table(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
    ncounts_file <- counts_file[,-1]
    final_data_new=ncounts_file[,grepl(paste0(manifest_file$samples, collapse = "|"), names(ncounts_file))]
    filt_data_new <- final_data_new[match(as.character(manifest_file$samples), names(final_data_new))]
    nfilt_data_new <- cbind(counts_file[,1], filt_data_new)
    if (input$geneid=="gene_id"){
      nfilt_data_new[,1] <- gsub('\\..*','',nfilt_data_new[,1])
      colnames(nfilt_data_new) <- c("gene_id", as.character(manifest_file$shortnames))
      housekeepers<- data.frame(gene_id= c("ENSG00000111640", "ENSG00000075624"))
      nfilt_data_new$gene_id <- toupper(nfilt_data_new$gene_id)
      housekeep <- merge(housekeepers, nfilt_data_new, by="gene_id", all=F)
      if (dim(housekeep)[1] == 2){
        barplot(as.matrix(log2(housekeep[,-1])), names.arg = as.character(manifest_file[,2]), las=2, col=c("orange", "blue"), ylim=c(0,50))
        legend('topright', legend=c("ENSG00000111640", "ENSG00000075624"), pch=15, col= c("orange", "blue"))
        
      }
      else {
        plot(NA, xlim=c(0,1), ylim=c(0,1))
        text(0.5,0.5,"No expression of ENSG00000111640 and ENSG00000075624 found", cex = 1)
      }	
    }
    else if (input$geneid=="gene_symbol"){
      colnames(nfilt_data_new) <- c("gene_symbol", as.character(manifest_file$shortnames))
      housekeepers<- data.frame(gene_symbol= c("GAPDH", "ACTB"))
      nfilt_data_new$gene_symbol <- toupper(nfilt_data_new$gene_symbol)
      housekeep <- merge(housekeepers, nfilt_data_new, by="gene_symbol", all=F)
      if (dim(housekeep)[1] == 2){
        barplot(as.matrix(log2(housekeep[,-1])), names.arg = as.character(manifest_file[,2]), las=2, col=c("orange", "blue"), ylim=c(0,50))
        legend('topright', legend=c("GAPDH", "ACTB"), pch=15, col= c("orange", "blue"))
      }
      else {
        plot(NA, xlim=c(0,1), ylim=c(0,1))
        text(0.5,0.5,"No expression of GAPDH and ACTB found", cex = 1)
    }
  }   
}


select_plot4 <- function(){
    manifest_file <- read.delim(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
    manifest_file <- manifest_file[manifest_file$include,]
    counts_file <- read.delim(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
    if (input$ftype=="Raw.Counts"){
    filt_data_new <- filtering_data(manifest_file, counts_file) 
    data.norm <- normalize_data(filt_data_new)
    data.final <- as.matrix(data.norm)
    colnames(data.final) <- manifest_file$shortnames
    respca <- PCA(data.final, graph = FALSE)
    p <- fviz_pca_var(respca, col.circle = rgb(255,255,255, alpha= 0, maxColorValue = 255), labelsize= 4, pointsize= 4, label= "var", geom = c("point", "text"), habillage=manifest_file$groups, title= "Principal Component Variances for all samples- zscored Normalized") + theme(text = element_text(size = 14))
    plot(p)
    }
    else if (input$ftype=="Normalized.Counts"){
      filt_data_new <- filtering_data(manifest_file, counts_file)
      data.norm <- normalize_data(filt_data_new)
      data.final <- as.matrix(data.norm)
      data.final <- filt_data_new
      colnames(data.final) <- manifest_file$shortnames
      respca <- PCA(data.final, graph = FALSE)
      p <- fviz_pca_var(respca, col.circle = rgb(255,255,255, alpha= 0, maxColorValue = 255), labelsize= 4, pointsize= 4, label= "var", geom = c("point", "text"), habillage=manifest_file$groups, title= "Principal Component Variances for all samples- zscored Normalized") + theme(text = element_text(size = 14))
      plot(p)
  }
  
}

select_plot5 <- function(){
    manifest_file <- read.delim(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
    manifest_file <- manifest_file[manifest_file$include,]
    counts_file <- read.delim(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
    if (input$ftype=="Raw.Counts"){
      filt_data_new <- filtering_data(manifest_file, counts_file) 
      data.final <- filt_data_new
      colnames(data.final) <- manifest_file$shortnames
      respca <- PCA(data.final, graph = FALSE)
      p <- fviz_pca_var(respca, col.circle = rgb(255,255,255, alpha= 0, maxColorValue = 255), labelsize= 4, pointsize= 4, label= "var", geom = c("point", "text"), habillage=manifest_file$groups, title= "Principal Component Variances for all samples- Un-Normalized") + theme(text = element_text(size = 14))
      plot(p)
    }
    else if (input$ftype=="Normalized.Counts"){
      filt_data_new <- filtering_data(manifest_file, counts_file) 
      data.final <- filt_data_new
      colnames(data.final) <- manifest_file$shortnames
      respca <- PCA(data.final, graph = FALSE)
      p <- fviz_pca_var(respca, col.circle = rgb(255,255,255, alpha= 0, maxColorValue = 255), labelsize= 4, pointsize= 4, label= "var", geom = c("point", "text"), habillage=manifest_file$groups, title= "Principal Component Variances for all samples- Un-Normalized") + theme(text = element_text(size = 14))
      plot(p)
  }
}

select_plot6 <- function(){
     manifest_file <- read.delim(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
    manifest_file <- manifest_file[manifest_file$include,]
    counts_file <- read.delim(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
    filt_data_new <- filtering_data(manifest_file, counts_file) 
	data.norm <- normalize_data(filt_data_new)
	data.final <- as.matrix(data.norm)
	librarysize <- colSums(data.final)
      if(dim(manifest_file)[2]==5){
        final.tran_data.final <- as.data.frame(cbind(as.numeric(manifest_file[,5]),librarysize))
        MFA_respca <- MFA(final.tran_data.final, group = c(1,1), type = c('s','s'), name.group = c(names(manifest_file[5]),"Lib.Size"), graph = F) 
        respca <- fviz_mfa_var(MFA_respca, "group", labelsize= 3, pointsize= 4, title= "Multifactor Principal Component Variances")+ theme(text = element_text(size = 14))
        plot(respca)
        
      }
      else if(dim(manifest_file)[2]==6){
        final.tran_data.final <- as.data.frame(cbind(as.numeric(manifest_file[,5]),as.numeric(manifest_file[,6]),librarysize))
        MFA_respca <- MFA(final.tran_data.final, group = c(1,1,1), type = c('s','s','s'), name.group = c(names(manifest_file[5]),names(manifest_file[6]),"Lib.Size"), graph = F) 
        respca <- fviz_mfa_var(MFA_respca, "group", labelsize= 4, pointsize= 4, title= "Multifactor Principal Component Variances")+ theme(text = element_text(size = 14))
        plot(respca)
      }
      else if(dim(manifest_file)[2]==7){
        final.tran_data.final <- as.data.frame(cbind(as.numeric(manifest_file[,5]),as.numeric(manifest_file[,6]),as.numeric(manifest_file[,7]),librarysize))
        MFA_respca <- MFA(final.tran_data.final, group = c(1,1,1,1), type = c('s','s','s','s'), name.group = c(names(manifest_file[5]),names(manifest_file[6]),names(manifest_file[7]),"Lib.Size"), graph = F) 
        respca <- fviz_mfa_var(MFA_respca, "group", labelsize= 4, pointsize= 4, title= "Multifactor Principal Component Variances")+ theme(text = element_text(size = 14))
        plot(respca)
      }
      else if(dim(manifest_file)[2]==8){
        final.tran_data.final <- as.data.frame(cbind(as.numeric(manifest_file[,5]),as.numeric(manifest_file[,6]),as.numeric(manifest_file[,7]),as.numeric(manifest_file[,8]),librarysize))
        MFA_respca <- MFA(final.tran_data.final, group = c(1,1,1,1,1), type = c('s','s','s','s','s'), name.group = c(names(manifest_file[5]),names(manifest_file[6]),names(manifest_file[7]),names(manifest_file[8]),"Lib.Size"), graph = F) 
        respca <- fviz_mfa_var(MFA_respca, "group", labelsize= 4, pointsize= 4, title= "Multifactor Principal Component Variances")+ theme(text = element_text(size = 14))
        plot(respca)
      }
      else if(dim(manifest_file)[2]==9){
        final.tran_data.final <- as.data.frame(cbind(as.numeric(manifest_file[,5]),as.numeric(manifest_file[,6]),as.numeric(manifest_file[,7]),as.numeric(manifest_file[,8]),as.numeric(manifest_file[,9]),librarysize))
        MFA_respca <- MFA(final.tran_data.final, group = c(1,1,1,1,1,1), type = c('s','s','s','s','s','s'), name.group = c(names(manifest_file[5]),names(manifest_file[6]),names(manifest_file[7]),names(manifest_file[8]),names(manifest_file[9]),"Lib.Size"), graph = F) 
        respca <- fviz_mfa_var(MFA_respca, "group", labelsize= 4, pointsize= 4, title= "Multifactor Principal Component Variances")+ theme(text = element_text(size = 14))
        plot(respca)
      }
      else if(dim(manifest_file)[2]==10){
        final.tran_data.final <- as.data.frame(cbind(as.numeric(manifest_file[,5]),as.numeric(manifest_file[,6]),as.numeric(manifest_file[,7]),as.numeric(manifest_file[,8]),as.numeric(manifest_file[,9]),as.numeric(manifest_file[,10]),librarysize))
        MFA_respca <- MFA(final.tran_data.final, group = c(1,1,1,1,1,1,1), type = c('s','s','s','s','s','s','s'), name.group = c(names(manifest_file[5]),names(manifest_file[6]),names(manifest_file[7]),names(manifest_file[8]),names(manifest_file[9]),names(manifest_file[10]),"Lib.Size"), graph = F) 
        respca <- fviz_mfa_var(MFA_respca, "group", labelsize= 4, pointsize= 4, title= "Multifactor Principal Component Variances")+ theme(text = element_text(size = 14))
        plot(respca)
      }
      else if(dim(manifest_file)[2]==11){
        final.tran_data.final <- as.data.frame(cbind(as.numeric(manifest_file[,5]),as.numeric(manifest_file[,6]),as.numeric(manifest_file[,7]),as.numeric(manifest_file[,8]),as.numeric(manifest_file[,9]),as.numeric(manifest_file[,10]),as.numeric(manifest_file[,11]),librarysize))
        MFA_respca <- MFA(final.tran_data.final, group = c(1,1,1,1,1,1,1,1), type = c('s','s','s','s','s','s','s','s'), name.group = c(names(manifest_file[5]),names(manifest_file[6]),names(manifest_file[7]),names(manifest_file[8]),names(manifest_file[9]),names(manifest_file[10]),names(manifest_file[11]),"Lib.Size"), graph = F) 
        respca <- fviz_mfa_var(MFA_respca, "group", labelsize= 4, pointsize= 4, title= "Multifactor Principal Component Variances")+ theme(text = element_text(size = 14))
        plot(respca)
      }
      else{
        plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
        text(0.5,0.5,"No factor to calculate", cex = 1.2)
  }
}

select_plot7 <- function(){
    manifest_file <- read.delim(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
    manifest_file <- manifest_file[manifest_file$include,]
    counts_file <- read.delim(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
    if (input$ftype=="Raw.Counts"){
    filt_data_new <- filtering_data(manifest_file, counts_file)
    filt_data_new <- normalize_data(filt_data_new)
    colnames(filt_data_new) <- manifest_file$shortnames
    phylo_plot(filt_data_new, 0.9)
    }
    else if (input$ftype=="Normalized.Counts"){
      filt_data_new <- filtering_data(manifest_file, counts_file)
      colnames(filt_data_new) <- manifest_file$shortnames
      phylo_plot(filt_data_new, 0.9)
  }
}

select_plot8 <- function(){
    manifest_file <- read.delim(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
    manifest_file <- manifest_file[manifest_file$include,]
    counts_file <- read.delim(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
    if (input$ftype=="Raw.Counts"){
    filt_data_new <- filtering_data(manifest_file, counts_file)
    cor_cts <- cor(filt_data_new, method="pearson")
    colnames(cor_cts) <- paste0(manifest_file$shortnames, ".", manifest_file$groups)
    row.names(cor_cts) <- paste0(manifest_file$shortnames, ".", manifest_file$groups)
    cols<- colorRampPalette(c("blue", "white", "red"))(20)
    corrplot(cor_cts, type="upper", order="hclust", col=cols, method = "square",tl.col="black", tl.srt=45, tl.cex = 1, main = "Pearson correlation of all the samples", mar=c(2,2,5,2))
    }
    else if (input$ftype=="Normalized.Counts"){
      filt_data_new <- filtering_data(manifest_file, counts_file)
      cor_cts <- cor(filt_data_new, method="pearson")
      colnames(cor_cts) <- paste0(manifest_file$shortnames, ".", manifest_file$groups)
      row.names(cor_cts) <- paste0(manifest_file$shortnames, ".", manifest_file$groups)
      cols<- colorRampPalette(c("blue", "white", "red"))(20)
      corrplot(cor_cts, type="upper", order="hclust", col=cols, method = "square",tl.col="black", tl.srt=45, tl.cex = 1, main = "Pearson correlation of all the samples", mar=c(2,2,5,2))
  }
  }

select_plot9 <- function(){
    manifest_file <- read.delim(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
    manifest_file <- manifest_file[manifest_file$include,]
    counts_file <- read.delim(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
    if (input$ftype=="Raw.Counts"){
    filt_data_new <- filtering_data(manifest_file, counts_file)
    cor_cts <- cor(filt_data_new, method="spearman")
    colnames(cor_cts) <- paste0(manifest_file$shortnames, ".", manifest_file$groups)
    row.names(cor_cts) <- paste0(manifest_file$shortnames, ".", manifest_file$groups)
    cols<- colorRampPalette(c("blue", "white", "red"))(20)
    corrplot(cor_cts, type="upper", order="hclust", col=cols, tl.col="black", tl.srt=45, tl.cex = 1,method = "square",main = "Spearman correlation of all the samples", mar=c(2,2,5,2))
    }
    else if (input$ftype=="Normalized.Counts"){
      filt_data_new <- filtering_data(manifest_file, counts_file)
      cor_cts <- cor(filt_data_new, method="spearman")
      colnames(cor_cts) <- paste0(manifest_file$shortnames, ".", manifest_file$groups)
      row.names(cor_cts) <- paste0(manifest_file$shortnames, ".", manifest_file$groups)
      cols<- colorRampPalette(c("blue", "white", "red"))(20)
      corrplot(cor_cts, type="upper", order="hclust", col=cols, tl.col="black", tl.srt=45, tl.cex = 1,method = "square",main = "Spearman correlation of all the samples", mar=c(2,2,5,2))
 }
  }

select_plot10 <- function(){
    manifest_file <- read.table(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
    manifest_file <- manifest_file[manifest_file$include,]
    counts_file <- read.table(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
    data_new = na.omit(counts_file)
    ndata_new= data_new[,-1]
    final_data_new=ndata_new[,grepl(paste0(manifest_file$samples, collapse = "|"), names(ndata_new))]
    filt_data_new <- final_data_new[match(as.character(manifest_file$samples), names(final_data_new))]
    data_new[,1] <- gsub('\\..*','',data_new[,1])
    row.names(filt_data_new) <- make.names(as.character(data_new[,1]), unique=T)
    colnames(filt_data_new) <- manifest_file$shortnames
    gene_query <- input$genes
    gene <- filt_data_new[grep(paste0("^",toupper(input$genes),"$"), toupper(row.names(filt_data_new))),]
    if (dim(gene)[1] == 1){
      barplot(as.matrix(log2(gene)), names.arg = as.character(manifest_file$shortnames), las=2, main= as.character(input$genes), ylab="Log2 expression", col="#fb6a4a", ylim=c(0,max(log2(gene[,-1]+5))))
    }
    else if (dim(gene)[1] == 0){
      plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
      text(0.5,0.5,"Please choose gene IDs", cex = 1.2)
    }
    else {
      plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
      text(0.5,0.5,"Cannot Plot- check if you are using right gene identifier ", cex = 1.2)
  }
}

select_plot11 <- function(){
    manifest_file <- read.table(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
    manifest_file <- manifest_file[manifest_file$include,]
    counts_file <- read.table(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
    ncts = na.omit(counts_file)
    ndata_new= ncts[,-1]
    final_data_new=ndata_new[,grepl(paste0(manifest_file$samples, collapse = "|"), names(ndata_new))]
    filt_data_new <- final_data_new[match(as.character(manifest_file$samples), names(final_data_new))]
    
    if(input$geneid=="gene_id"){
      row.names(filt_data_new) <- make.names(gsub('\\..*','',ncts[,1]), unique=T)
      log.filt_data_new <- log(filt_data_new)
      filt.log.filt_data_new <- log.filt_data_new[!is.infinite(rowSums(log.filt_data_new)),]
      filt.log.filt_data_new =setDT(filt.log.filt_data_new, keep.rownames = T)[]
      names(filt.log.filt_data_new)[1]<- "gene_id"
      
      if(input$org=="Human"){
        nhg <- read.table('datafiles/annotation_hg38withGC.txt', sep = '\t', header = T, check.names = F)
      }
      else if(input$org=="Mouse"){
        nhg <- read.table('datafiles/annotation_m38withGC.txt', sep = '\t', header = T, check.names = F)
      }
      
      else{
        nhg <- data.frame(a=as.Date(character()),
                          b=character(), 
                          c=character(), 
                          stringsAsFactors=FALSE) 
      }
      
      if (dim(nhg)[1]>1){
        filt_data_new_withgc <- merge(nhg, filt.log.filt_data_new, by="gene_id")
        final_len <- as.matrix(filt_data_new_withgc$length)
        final_gc <- as.matrix(filt_data_new_withgc$gc)
        final_data <- as.matrix(filt_data_new_withgc[,-c(1,2,3,4)])
        cl <- rainbow(dim(final_data)[2])
        plot(NA, xlim=c(0.1,0.9),ylim=c(1,10), xlab="GC content", ylab="Gene counts")
        plotcol=NA
        for (i in 1:dim(final_data)[2]){
          points(loess.smooth(final_gc,final_data[,i]), type='l', col=cl[i],lwd=1.2)
          plotcol[i] <- cl[i]
        }
        legend("top", legend = manifest_file$shortnames, col = plotcol, lwd = 1.2, cex = 1, bty='n', xpd=T, ncol=5)
      }
      else {
        plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
        text(0.5,0.5,"Cannot Plot- the Organism you provided is not compatible", cex = 1.2)
      }
      }
      

      else if(input$geneid=="gene_symbol"){
        row.names(filt_data_new) <- make.names(ncts[,1], unique=T)
        log.filt_data_new <- log(filt_data_new)
        filt.log.filt_data_new <- log.filt_data_new[!is.infinite(rowSums(log.filt_data_new)),]
        filt.log.filt_data_new =setDT(filt.log.filt_data_new, keep.rownames = T)[]
        names(filt.log.filt_data_new)[1]<- "gene_symbol"
        
        if(input$org=="Human"){
          nhg <- read.table('datafiles/annotation_hg38withGC.txt', sep = '\t', header = T, check.names = F)
        }
        else if(input$org=="Mouse"){
          nhg <- read.table('datafiles/annotation_m38withGC.txt', sep = '\t', header = T, check.names = F)
        }
        else{
          nhg <- data.frame(a=as.Date(character()),
                            b=character(), 
                            c=character(), 
                            stringsAsFactors=FALSE) 
        }
        
        if (dim(nhg)[1]>1){
          
          filt_data_new_withgc <- merge(nhg, filt.log.filt_data_new, by="gene_symbol")
          final_len <- as.matrix(filt_data_new_withgc$length)
          final_gc <- as.matrix(filt_data_new_withgc$gc)
          final_data <- as.matrix(filt_data_new_withgc[,-c(1,2,3,4)])
          cl <- rainbow(dim(final_data)[2])
          plot(NA, xlim=c(0.1,0.9),ylim=c(1,10), xlab="GC content", ylab="Gene counts")
          plotcol=NA
          for (i in 1:dim(final_data)[2]){
            points(loess.smooth(final_gc,final_data[,i]), type='l', col=cl[i],lwd=1.2)
            plotcol[i] <- cl[i]
          }
          legend("top", legend = manifest_file$shortnames, col = plotcol, lwd = 1.2, cex = 1, bty='n', xpd=T, ncol=5)
        }
        
        else {
          plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',xaxt='n',yaxt='n', xlab = NA, ylab = NA)
          text(0.5,0.5,"Cannot Plot- the Organism you provided is not compatible", cex = 1.2)
        }

      }
}


output$table1.output <- renderPrint({ sum_input()})
output$plot_sum2.output <- renderPlot({ sum2_input()})
output$plot_sum3.output <- renderPlot({ sum3_input()})
output$plot_sum4.output <- renderPlot({ sum4_input()})
output$pca2.output <- renderPlot({pca2_input()})
output$pca3.output <- renderPlot({pca3_input()})
output$pca4.output <- renderPlot({pca4_input()})
output$cor1.output <- renderPlot({cor1_input()})
output$cor2.output <- renderPlot({cor2_input()})
output$cor3.output <- renderPlot({cor3_input()})
output$exp.output <- renderPlot({exp_input()})
output$bias.output <- renderPlot({bias_input()})


output$sum2download.output <- downloadHandler(filename= function(){"countsDistribution_plot.pdf"}, content= function(file){pdf(file);select_plot1();dev.off()})
output$sum3download.output <- downloadHandler(filename= function(){"detectedtranscripts_plot.pdf"}, content= function(file){pdf(file);select_plot2();dev.off()})
output$sum4download.output <- downloadHandler(filename= function(){"housekeeping_plot.pdf"}, content= function(file){pdf(file);select_plot3();dev.off()})
output$pca2download.output <- downloadHandler(filename= function(){"PCA_plot_normalized.pdf"}, content= function(file){pdf(file);select_plot4();dev.off()})
output$pca3download.output <- downloadHandler(filename= function(){"PCA_plot_non-normalized.pdf"}, content= function(file){pdf(file);select_plot5();dev.off()})
output$pca4download.output <- downloadHandler(filename= function(){"PCAmultifactor_plot.pdf"}, content= function(file){pdf(file);select_plot6();dev.off()})
output$cor1download.output <- downloadHandler(filename= function(){"Hierarchical_plot.pdf"}, content= function(file){pdf(file);select_plot7();dev.off()})
output$cor2download.output <- downloadHandler(filename= function(){"PearsonCorrelation_plot.pdf"}, content= function(file){pdf(file);select_plot8();dev.off()})
output$cor3download.output <- downloadHandler(filename= function(){"SpearmanCorrelation_plot.pdf"}, content= function(file){pdf(file);select_plot9();dev.off()})
output$expdownload.output <- downloadHandler(filename= function(){"Expression_plot.pdf"}, content= function(file){pdf(file);select_plot10();dev.off()})
output$biasdownload.output <- downloadHandler(filename= function(){"GCbias_plot.pdf"}, content= function(file){pdf(file);select_plot11();dev.off()})
output$downloadmanifest <- downloadHandler(filename= function(){'samplemanifestfile.txt'}, content= function(file){file.copy("datafiles/samplemanifestfile.txt",file)})
output$downloadmanifestwithmfa <- downloadHandler(filename= function(){'samplemanifestfile_multifactor.txt'}, content= function(file){file.copy("datafiles/samplemanifestfile_multifactor.txt",file)})
output$downloadrawcountsgeneid <- downloadHandler(filename= function(){'geneid_rawcounts.txt'}, content= function(file){file.copy("datafiles/geneid_rawcounts.txt",file)})
output$downloadrawcountsgenesymbol <- downloadHandler(filename= function(){'genesymbol_rawcounts.txt'}, content= function(file){file.copy("datafiles/genesymbol_rawcounts.txt",file)})

})