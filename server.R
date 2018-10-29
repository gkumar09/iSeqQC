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
source("R/correlation.R")
source("R/density.R")
source("R/summary.R")

shinyServer(function(input, output, session) {


sum_input <- reactive({
    cts_file <- input$countsmatrix
    if (is.null(cts_file)){
      return(invisible())}
    else{
      manifest_file <- read.delim(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
      manifest_file <- manifest_file[manifest_file$include,]
      counts_file <- read.delim(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
      filt_data_new <- filtering_data(manifest_file, counts_file)
      colnames(filt_data_new) <- manifest_file$shortnames
      data_final <- SumMarize(filt_data_new)
      data_final
    }
})

sum2_input <- reactive({
  cts_file <- input$countsmatrix
  if (is.null(cts_file)){
    return(invisible())}
  else{
    manifest_file <- read.delim(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
    manifest_file <- manifest_file[manifest_file$include,]
    counts_file <- read.delim(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
    filt_data_new <- filtering_data(manifest_file, counts_file)
    colnames(filt_data_new) <- manifest_file$shortnames
    distributionplot(filt_data_new)
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
    filt_data_new <- filtering_data(manifest_file, counts_file)
    colnames(filt_data_new) <- manifest_file$shortnames
    transcriptplot(filt_data_new)
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
    filt_data_new <- filtering_data(manifest_file, counts_file) 
    data.norm <- normalize_data(filt_data_new)
    data.final <- as.matrix(data.norm)
    colnames(data.final) <- manifest_file$shortnames
    respca <- PCA(data.final, graph = FALSE)
    p <- fviz_pca_var(respca, col.circle = rgb(255,255,255, alpha= 0, maxColorValue = 255), labelsize= 4, pointsize= 4, label= "var", geom = c("point", "text"), habillage=manifest_file$groups, title= "Principal Component Variances for all samples") + theme(text = element_text(size = 14))
    plot(p)
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
    filt_data_new <- filtering_data(manifest_file, counts_file)
    data.norm <- normalize_data(filt_data_new)
    colnames(data.norm) <- manifest_file$shortnames
    phylo_plot(data.norm, 0.9)
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
    filt_data_new <- filtering_data(manifest_file, counts_file)
    cor_cts <- cor(filt_data_new, method="pearson")
    colnames(cor_cts) <- paste0(manifest_file$shortnames, ".", manifest_file$groups)
    row.names(cor_cts) <- paste0(manifest_file$shortnames, ".", manifest_file$groups)
    cols<- colorRampPalette(c("blue", "white", "red"))(20)
    corrplot(cor_cts, type="upper", order="hclust", col=cols, addCoef.col = "black", method = "color", tl.col="black", tl.srt=45, tl.cex = 0.6, main = "Pearson correlation of all the samples", mar=c(2,2,5,2))
  }
})


# downloading plots
select_plot1 <- function(){
  manifest_file <- read.delim(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
  manifest_file <- manifest_file[manifest_file$include,]
  counts_file <- read.delim(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
  filt_data_new <- filtering_data(manifest_file, counts_file)
  colnames(filt_data_new) <- manifest_file$shortnames
  distributionplot(filt_data_new)
}


select_plot2 <- function(){
  manifest_file <- read.delim(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
  manifest_file <- manifest_file[manifest_file$include,]
  counts_file <- read.delim(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
  filt_data_new <- filtering_data(manifest_file, counts_file)
  colnames(filt_data_new) <- manifest_file$shortnames
  transcriptplot(filt_data_new)
}


select_plot4 <- function(){
  manifest_file <- read.delim(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
  manifest_file <- manifest_file[manifest_file$include,]
  counts_file <- read.delim(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
  filt_data_new <- filtering_data(manifest_file, counts_file) 
  data.norm <- normalize_data(filt_data_new)
  data.final <- as.matrix(data.norm)
  colnames(data.final) <- manifest_file$shortnames
  respca <- PCA(data.final, graph = FALSE)
  p <- fviz_pca_var(respca, col.circle = rgb(255,255,255, alpha= 0, maxColorValue = 255), labelsize= 4, pointsize= 4, label= "var", geom = c("point", "text"), habillage=manifest_file$groups, title= "Principal Component Variances for all samples") + theme(text = element_text(size = 14))
  plot(p)
}

select_plot6 <- function(){
  manifest_file <- read.delim(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
  manifest_file <- manifest_file[manifest_file$include,]
  counts_file <- read.delim(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
  filt_data_new <- filtering_data(manifest_file, counts_file)
  data.norm <- normalize_data(filt_data_new)
  colnames(data.norm) <- manifest_file$shortnames
  phylo_plot(data.norm, 0.4)
}

select_plot7 <- function(){
  manifest_file <- read.delim(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
  manifest_file <- manifest_file[manifest_file$include,]
  counts_file <- read.delim(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
  filt_data_new <- filtering_data(manifest_file, counts_file)
  cor_cts <- cor(filt_data_new, method="pearson")
  print(cor_cts)
  colnames(cor_cts) <- paste0(manifest_file$shortnames, ".", manifest_file$groups)
  row.names(cor_cts) <- paste0(manifest_file$shortnames, ".", manifest_file$groups)
  cols<- colorRampPalette(c("blue", "white", "red"))(20)
  corrplot(cor_cts, type="upper", order="hclust", col=cols, addCoef.col = "black", method = "color", tl.col="black", tl.srt=45, tl.cex = 0.6, main = "Pearson correlation of all the samples", mar=c(2,2,5,2))
}


output$table1.output <- renderPrint({ sum_input()})
output$plot_sum2.output <- renderPlot({ sum2_input()})
output$plot_sum3.output <- renderPlot({ sum3_input()})
output$plot2.output <- renderPlot({pca2_input()})
output$plot4.output <- renderPlot({cor1_input()})
output$plot5.output <- renderPlot({cor2_input()})


output$sum1download.output <- downloadHandler(filename= function(){"countsDistribution_plot.pdf"}, content= function(file){pdf(file);select_plot1();dev.off()})
output$sum2download.output <- downloadHandler(filename= function(){"detectedtranscripts_plot.pdf"}, content= function(file){pdf(file);select_plot2();dev.off()})
output$pca2download.output <- downloadHandler(filename= function(){"PCA_plot.pdf"}, content= function(file){pdf(file);select_plot4();dev.off()})
output$hierdownload.output <- downloadHandler(filename= function(){"Hierarchical_plot.pdf"}, content= function(file){pdf(file);select_plot6();dev.off()})
output$cordownload.output <- downloadHandler(filename= function(){"Correlational_plot.pdf"}, content= function(file){pdf(file);select_plot7();dev.off()})
})