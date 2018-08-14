options(shiny.maxRequestSize = 1000000*1024^2)

source("functions/correlation.R")
source("functions/density.R")
source("functions/distribution.R")
source("functions/filtering.R")
source("functions/housekeeping.R")
source("functions/housekeeping_filt.R")
source("functions/mds.R")
source("functions/normalizing.R")
source("functions/pca.R")
source("functions/phylo.R")
source("functions/summary.R")
source("functions/weight.R")

library("shiny")
library("FactoMineR")
library("factoextra")
library("som")
library("psych")
library("data.table")
library("ape")
library("corrplot")
library("limma")

shinyServer(function(input, output, session) {


sum_input <- reactive({
    cts_file <- input$countsmatrix
    if (is.null(cts_file)){
      return(invisible())}
    else{
      manifest_file <- read.table(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
      manifest_file <- manifest_file[manifest_file$include,]
      if (any(names(manifest_file) != c("samples","shortnames","groups", "include"))){
        err= "FALSE"
        validate(err != "FALSE", "Sorry, unable to plot: Either the column names in sample phenotype file doesn't meet the recommended names or order")
      }
      counts_file <- read.table(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
      if (any(grepl(paste0(manifest_file$samples, collapse = "|"), names(counts_file))[-1] == "FALSE")){
        err= "FALSE"
        validate(err != "FALSE", "Sorry, can't plot: The sample names in counts file doesn't match with 'samples' in sample phenotype file")
      }
      else {
      filt_data_new <- filtering_data(manifest_file, counts_file)
      colnames(filt_data_new) <- manifest_file$shortnames
      data_final <- SumMarize(filt_data_new)
      data_final
      }
    }
})

sum2_input <- reactive({
  cts_file <- input$countsmatrix
  if (is.null(cts_file)){
    return(invisible())}
  else{
    manifest_file <- read.table(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
    manifest_file <- manifest_file[manifest_file$include,]
    if (any(names(manifest_file) != c("samples","shortnames","groups", "include"))){
      err= "FALSE"
      validate(err != "FALSE", "Sorry, unable to plot: Either the column names in sample phenotype file doesn't meet the recommended names or order")
    }
    counts_file <- read.table(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
    if (any(grepl(paste0(manifest_file$samples, collapse = "|"), names(counts_file))[-1] == "FALSE")){
      err= "FALSE"
      validate(err != "FALSE", "Sorry, unable to plot: The sample names in counts file doesn't match with 'samples' in sample phenotype file")
    }
    else {
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
    manifest_file <- read.table(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
    manifest_file <- manifest_file[manifest_file$include,]
    if (any(names(manifest_file) != c("samples","shortnames","groups", "include"))){
      err= "FALSE"
      validate(err != "FALSE", "Sorry, unable to plot: Either the column names in sample phenotype file doesn't meet the recommended names or order")
    }
    counts_file <- read.table(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
    if (any(grepl(paste0(manifest_file$samples, collapse = "|"), names(counts_file))[-1] == "FALSE")){
      err= "FALSE"
      validate(err != "FALSE", "Sorry, unable to plot: The sample names in counts file doesn't match with 'samples' in sample phenotype file")
    }
    else {
    filt_data_new <- filtering_data(manifest_file, counts_file)
    colnames(filt_data_new) <- manifest_file$shortnames
    transcriptplot(filt_data_new)
    }
  }
})


pca1_input <- reactive({
    cts_file <- input$countsmatrix
    if (is.null(cts_file)){
      return(NULL)}
    else{
      manifest_file <- read.table(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
      manifest_file <- manifest_file[manifest_file$include,]
      if (any(names(manifest_file) != c("samples","shortnames","groups", "include"))){
        err= "FALSE"
        validate(err != "FALSE", "Sorry, unable to plot: Either the column names in sample phenotype file doesn't meet the recommended names or order")
      }
      counts_file <- read.table(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
      if (any(grepl(paste0(manifest_file$samples, collapse = "|"), names(counts_file))[-1] == "FALSE")){
        err= "FALSE"
        validate(err != "FALSE", "Sorry, unable to plot: The sample names in counts file doesn't match with 'samples' in sample phenotype file")
      }
      else {
      filt_data_new <- filtering_data(manifest_file, counts_file) 
      data.norm <- normalize_data(filt_data_new)
      mdsplot(data.norm, manifest_file, 1.0)
      }
    }
})

pca2_input <- reactive({
  cts_file <- input$countsmatrix
  if (is.null(cts_file)){
    return(NULL)}
  else{
    manifest_file <- read.table(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
    manifest_file <- manifest_file[manifest_file$include,]
    if (any(names(manifest_file) != c("samples","shortnames","groups", "include"))){
      err= "FALSE"
      validate(err != "FALSE", "Sorry, unable to plot: Either the column names in sample phenotype file doesn't meet the recommended names or order")
    }
    counts_file <- read.table(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
    if (any(grepl(paste0(manifest_file$samples, collapse = "|"), names(counts_file))[-1] == "FALSE")){
      err= "FALSE"
      validate(err != "FALSE", "Sorry, unable to plot: The sample names in counts file doesn't match with 'samples' in sample phenotype file")
    }
    else {
    filt_data_new <- filtering_data(manifest_file, counts_file) 
    data.norm <- normalize_data(filt_data_new)
    data.final <- as.matrix(data.norm)
    colnames(data.final) <- manifest_file$shortnames
    pcaplot(data.final, manifest_file)
    }
  }
})

pca3_input <- reactive({
  cts_file <- input$countsmatrix
  if (is.null(cts_file)){
    return(NULL)}
  else{
    manifest_file <- read.table(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
    manifest_file <- manifest_file[manifest_file$include,]
    if (any(names(manifest_file) != c("samples","shortnames","groups", "include"))){
      err= "FALSE"
      validate(err != "FALSE", "Sorry, unable to plot: Either the column names in sample phenotype file doesn't meet the recommended names or order")
    }
    manifest_file <- manifest_file[order(manifest_file$groups),]
    counts_file <- read.table(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
    if (any(grepl(paste0(manifest_file$samples, collapse = "|"), names(counts_file))[-1] == "FALSE")){
      err= "FALSE"
      validate(err != "FALSE", "Sorry, unable to plot: The sample names in counts file doesn't match with 'samples' in sample phenotype file")
    }
    else {
    data = na.omit(counts_file)
    new_counts <- data[,-1]
    norma <- normalize_data(new_counts)
    new_norm <- data.frame(cbind(as.character(data[,1]), norma))
    names(new_norm)<- as.character(names(data))
    housekeeper <- housekeeping_filt(counts_file, new_norm)
    filt_data_new <- filtering_data(manifest_file, housekeeper)
    colnames(filt_data_new) <- paste0(manifest_file$shortnames, "_", manifest_file$groups)
    housekeeping_plot(manifest_file, filt_data_new, housekeeper, 0.7)    
    }
  }
})

cor1_input <- reactive({
  cts_file <- input$countsmatrix
  if (is.null(cts_file)){
    return(NULL)}
  else{
    manifest_file <- read.table(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
    manifest_file <- manifest_file[manifest_file$include,]
    if (any(names(manifest_file) != c("samples","shortnames","groups", "include"))){
      err= "FALSE"
      validate(err != "FALSE", "Sorry, unable to plot: Either the column names in sample phenotype file doesn't meet the recommended names or order")
    }
    counts_file <- read.table(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
    if (any(grepl(paste0(manifest_file$samples, collapse = "|"), names(counts_file))[-1] == "FALSE")){
      err= "FALSE"
      validate(err != "FALSE", "Sorry, unable to plot: The sample names in counts file doesn't match with 'samples' in sample phenotype file")
    }
    else {
    filt_data_new <- filtering_data(manifest_file, counts_file)
    data.norm <- normalize_data(filt_data_new)
    colnames(data.norm) <- paste0(manifest_file$shortnames, "_", manifest_file$groups)
    phylo_plot(data.norm, 0.9)
    }
  }
})

cor2_input <- reactive({
  cts_file <- input$countsmatrix
  if (is.null(cts_file)){
    return(NULL)}
  else{
    manifest_file <- read.table(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
    manifest_file <- manifest_file[manifest_file$include,]
    if (any(names(manifest_file) != c("samples","shortnames","groups", "include"))){
      err= "FALSE"
      validate(err != "FALSE", "Sorry, unable to plot: Either the column names in sample phenotype file doesn't meet the recommended names or order")
    }
    counts_file <- read.table(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
    if (any(grepl(paste0(manifest_file$samples, collapse = "|"), names(counts_file))[-1] == "FALSE")){
      err= "FALSE"
      validate(err != "FALSE", "Sorry, unable to plot: The sample names in counts file doesn't match with 'samples' in sample phenotype file")
    }
    else {
    filt_data_new <- filtering_data(manifest_file, counts_file)
    data.norm <- normalize_data(filt_data_new)
    colnames(data.norm) <- paste0(manifest_file$shortnames, ".", manifest_file$groups)
    correlation_plot(data.norm)
    }
  }
})

cor3_input <- reactive({
  cts_file <- input$countsmatrix
  if (is.null(cts_file)){
    return(NULL)}
  else{
    manifest_file <- read.table(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
    manifest_file <- manifest_file[manifest_file$include,]
    if (any(names(manifest_file) != c("samples","shortnames","groups", "include"))){
      err= "FALSE"
      validate(err != "FALSE", "Sorry, unable to plot: Either the column names in sample phenotype file doesn't meet the recommended names or order")
    }
    counts_file <- read.table(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
    if (any(grepl(paste0(manifest_file$samples, collapse = "|"), names(counts_file))[-1] == "FALSE")){
      err= "FALSE"
      validate(err != "FALSE", "Sorry, unable to plot: The sample names in counts file doesn't match with 'samples' in sample phenotype file")
    }
    else {
    filt_data_new <- filtering_data(manifest_file, counts_file)
    colnames(filt_data_new) <- paste0(manifest_file$shortnames, "_", manifest_file$groups)
    weight_plot(manifest_file, filt_data_new)
    }
  }
})

# downloading plots
select_plot1 <- function(){
  manifest_file <- read.table(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
  manifest_file <- manifest_file[manifest_file$include,]
  counts_file <- read.table(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
  filt_data_new <- filtering_data(manifest_file, counts_file)
  colnames(filt_data_new) <- manifest_file$shortnames
  distributionplot(filt_data_new)
}


select_plot2 <- function(){
  manifest_file <- read.table(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
  manifest_file <- manifest_file[manifest_file$include,]
  counts_file <- read.table(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
  filt_data_new <- filtering_data(manifest_file, counts_file)
  colnames(filt_data_new) <- manifest_file$shortnames
  transcriptplot(filt_data_new)
}

select_plot3 <- function(){
  manifest_file <- read.table(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
  manifest_file <- manifest_file[manifest_file$include,]
  counts_file <- read.table(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
  filt_data_new <- filtering_data(manifest_file, counts_file) 
  data.norm <- normalize_data(filt_data_new)
  mdsplot(data.norm, manifest_file, 0.6)
}

select_plot4 <- function(){
  manifest_file <- read.table(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
  manifest_file <- manifest_file[manifest_file$include,]
  counts_file <- read.table(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
  filt_data_new <- filtering_data(manifest_file, counts_file) 
  data.norm <- normalize_data(filt_data_new)
  data.final <- as.matrix(data.norm)
  colnames(data.final) <- manifest_file$shortnames
  pca <- pcaplot(data.final, manifest_file)
  plot(pca)
}

select_plot5 <- function(){
  manifest_file <- read.table(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
  manifest_file <- manifest_file[manifest_file$include,]
  manifest_file <- manifest_file[order(manifest_file$groups),]
  counts_file <- read.table(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
  data = na.omit(counts_file)
  new_counts <- data[,-1]
  norma <- normalize_data(new_counts)
  new_norm <- data.frame(cbind(as.character(data[,1]), norma))
  names(new_norm)<- as.character(names(data))
  housekeeper <- housekeeping_filt(counts_file, new_norm)
  filt_data_new <- filtering_data(manifest_file, housekeeper)
  colnames(filt_data_new) <- paste0(manifest_file$shortnames, "_", manifest_file$groups)
  par(mar=c(3,4,2,4))
  housekeeping_plot(manifest_file, filt_data_new, housekeeper, 0.4)
}

select_plot6 <- function(){
  manifest_file <- read.table(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
  manifest_file <- manifest_file[manifest_file$include,]
  counts_file <- read.table(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
  filt_data_new <- filtering_data(manifest_file, counts_file)
  data.norm <- normalize_data(filt_data_new)
  colnames(data.norm) <- paste0(manifest_file$shortnames, ".", manifest_file$groups)
  phylo_plot(data.norm, 0.6)
}

select_plot7 <- function(){
  manifest_file <- read.table(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
  manifest_file <- manifest_file[manifest_file$include,]
  counts_file <- read.table(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
  filt_data_new <- filtering_data(manifest_file, counts_file)
  data.norm <- normalize_data(filt_data_new)
  colnames(data.norm) <- paste0(manifest_file$shortnames, ".", manifest_file$groups)
  correlation_plot(data.norm)
}

select_plot8 <- function(){
  manifest_file <- read.table(input$manifest$datapath, sep="\t", header=TRUE, check.names = F)
  manifest_file <- manifest_file[manifest_file$include,]
  counts_file <- read.table(input$countsmatrix$datapath, sep="\t", header=TRUE, check.names = F)
  filt_data_new <- filtering_data(manifest_file, counts_file)
  colnames(filt_data_new) <- paste0(manifest_file$shortnames, "_", manifest_file$groups)
  weight_plot(manifest_file, filt_data_new)
}

output$table1.output <- renderPrint({ sum_input()})
output$plot_sum2.output <- renderPlot({ sum2_input()})
output$plot_sum3.output <- renderPlot({ sum3_input()})
output$plot1.output <- renderPlot({ Sys.sleep(2); pca1_input()})
output$plot2.output <- renderPlot({pca2_input()})
output$plot3.output <- renderPlot({pca3_input()})
output$plot4.output <- renderPlot({cor1_input()})
output$plot5.output <- renderPlot({cor2_input()})
output$plot6.output <- renderPlot({cor3_input()})

output$sum1download.output <- downloadHandler(filename= function(){"countsDistribution_plot.pdf"}, content= function(file){pdf(file);select_plot1();dev.off()})
output$sum2download.output <- downloadHandler(filename= function(){"detectedtranscripts_plot.pdf"}, content= function(file){pdf(file);select_plot2();dev.off()})
output$pca1download.output <- downloadHandler(filename= function(){"mds_plot.pdf"}, content= function(file){pdf(file);select_plot3();dev.off()})
output$pca2download.output <- downloadHandler(filename= function(){"PCA_plot.pdf"}, content= function(file){pdf(file);select_plot4();dev.off()})
output$pca3download.output <- downloadHandler(filename= function(){"housekeeping_average_plot.pdf"}, content= function(file){pdf(file);select_plot5();dev.off()})
output$hierdownload.output <- downloadHandler(filename= function(){"Hierarchical_plot.pdf"}, content= function(file){pdf(file);select_plot6();dev.off()})
output$cordownload.output <- downloadHandler(filename= function(){"Correlational_plot.pdf"}, content= function(file){pdf(file);select_plot7();dev.off()})
output$wtdownload.output <- downloadHandler(filename= function(){"sampleWeights_plot.pdf"}, content= function(file){pdf(file);select_plot8();dev.off()})
})