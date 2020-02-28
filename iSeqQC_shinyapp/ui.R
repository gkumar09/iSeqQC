library("shiny")


shinyUI(fluidPage
        (theme = "bootstrap.css",
        titlePanel(h2("iSeqQC- An Expression based Quality Control tool", style= "font-family: 'American Typewriter'; color:#081d58"), windowTitle = "iSeqQC"),
        hr(),
          navbarPage("",
                     tabPanel(h4("About"), 
                              titlePanel(h3("Welcome to iSeqQC",align="center",style = "font-family: 'times'; color:#081d58; font-si16pt; line-height:1.8", br(), 
                                            h4("It is an expression-based quality control tool to detect outliers either produced by batch effects or merely due to dissimilarity within a phenotypic group.", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.8"),br(),
                                            
                                            h4("Input: 1) Choose the file-type between raw or normalized counts options; ", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.8"),
                                            h4("2) Choose the organism of study;", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.8"),
                                            h4("3) Choose between gene symbol or gene id as in your first column of counts file", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.8"),
                                            h4("4) Sample phenotype file in tab-delimited with a minimum of 4 columns, 1) sample names; 2) short name of samples; 3) groups; 4) samples to include in the analysis; 5) any factors that could influence any bias such as library protocol or preparation method. The first 4 columns are mandatory and should strictly match the names and order as mentioned below (names case-sensitive), column 5-11 are optional and could be left blank", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.8"),br(),
                                            img(src='samplemanifest.png', width='400'),
                                            h4("5) A counts matrix file in tab-delimited file. This could be a raw or normalized (TPM/FPKM/log2/RPKM/RPM etc.). This file should have gene identifier (either gene symbol or gene id's) in first column as shown below", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.8"),br(),
                                            img(src='gs_cm.png', width='600'),br(), br(),
                                            h4("or",style = "font-family: 'times'; font-si14pt; line-height:1.4"),
                                            img(src='gi_cm.png', width='600'), br(),br(),br(),
                                            
                                            h4("Output: It is presented in two tabs, where in ‘Quality Control’, you are required to input the above-mentioned information and obtain output is in form of an expression summary table and several QC plots such as density plots, PCA, hierchichal clustering, correlation plots. In ‘Expression Plots’ tab, you can choose your gene of interest to plot expression profile among different phenotypic groups. ", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.8"),br(),
                                            
                                            h4("Please note: To remove any sample from analysis, simply replace 'TRUE' with 'FALSE' in 'include' column in sample manifest file.", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.4"), br(),
                                            h4("Additional information can be obtained at https://github.com/gkumar09/iSeqQC", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.8"),br(),
                                            
                                            h4("Citation:", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.4"),
                                            h4("Kumar G, Ertel A, Feldman G, Kupper J, Fortina P (2020). iSeqQC: A Tool for Expression-Based Quality Control in RNA Sequencing. BMC Bioinformatics. Feb 13;21(1):56. doi: 10.1186/s12859-020-3399-8. PMID: 32054449; PMCID: PMC7020508", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.2"),br(),
                                            h4("Contact:", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.4"),
                                            h4("Gaurav Kumar, PhD", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.2"),
                                            h4("Email address: gaurav.kumar@jefferson.edu", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.2"),
                                            br(),
                                            br()
                                             ))),
            tabPanel(h4("Quality Control"),  
                     sidebarPanel(selectInput('ftype', 'File-type', choices = c("counts-type","Raw.Counts","Normalized.Counts"), selectize=FALSE), selectInput('org', 'Organism', choices = c("Organism","Human","Mouse","Others"), selectize=FALSE), selectInput('geneid',"Gene Column identifier" ,choices = c("Gene-Identifier","gene_id","gene_symbol"), selectize=FALSE)),
                     sidebarPanel(
                     fileInput("manifest", label="Choose sample phenotype File", accept= "text/plain"), helpText("Please remember to use correct column names and matched sample names to avoid any errors")),
                     sidebarPanel(
                      tags$head(tags$style(type="text/css", "#loadmessage {position: fixed;top: 350px;left: 0px;width: 100%;padding: 5px 0px 5px 0px;text-align: center;font-weight: bold;font-size: 100%;color: #000000; z-index: 105;}"),
                               conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                               tags$div("Please wait while the server is working...",id="loadmessage"))),
                     tags$head(
                       tags$style(HTML(" .shiny-output-error-validation {color: red;}"))
                       ),
                     fileInput("countsmatrix", label="Choose Counts matrix File", accept= "text/plain"), helpText("Please remember to have official gene symbols or IDs to avoid any errors")),
                     #sidebarPanel(selectInput('org', 'Organism', choices = c("Human","Mouse","Others"), selectize=FALSE)),
                     
                     
                     column(12, align="center", br(), verbatimTextOutput(outputId= 'table1.output'), br(), br(), plotOutput('plot_sum2.output', width = "40%"),downloadButton(outputId = "sum2download.output", label = "Download the distribution of counts per sample"), br(), br(), plotOutput('plot_sum3.output', width = "50%"),downloadButton(outputId = "sum3download.output", label = "Download the Mapped Reads density per sample"),br(), br(),plotOutput('plot_sum4.output', width = "50%"),downloadButton(outputId = "sum4download.output", label = "Download the housekeeping genes expression per sample"), br(),br(),br(), br(),
                            plotOutput('pca2.output', width = "80%"),downloadButton(outputId = "pca2download.output", label = "Download the Principal components plot for all samples- normalized"), br(), br(), br(),plotOutput('pca3.output', width = "80%"),downloadButton(outputId = "pca3download.output", label = "Download the Principal components plot for all samples- non-normalized"),br(),br(), plotOutput('pca4.output', width = "80%"),downloadButton(outputId = "pca4download.output", label = "Download the Multifactor Principal components plot"),br(),br(),
                            plotOutput('cor1.output', width = "80%", height = "400px"),downloadButton(outputId = "cor1download.output", label = "Download the hierarchichal plot"), br(), br(), plotOutput('cor2.output', width = "80%", height = "500px"),downloadButton(outputId = "cor2download.output", label = "Download the pearson correlation plot"), br(), br(),plotOutput('cor3.output', width = "80%", height = "500px"),downloadButton(outputId = "cor3download.output", label = "Download the spearman correlation plot"),br(), br(),
                            plotOutput('bias.output', width = "80%", height = "400px"),downloadButton(outputId = "biasdownload.output", label = "Download the GC bias plot"), br(), br(),br()
                     )),
            tabPanel(h4("Expression Plots"),
                       
                     sidebarPanel(textInput("genes",label="Gene Name or ID"),helpText("Example: Select gene_id- if your file first column looks like ENSG00000111640Select gene_symbol- if your file first column looks like GAPDH")),
                     column(12, align="center", br(), plotOutput('exp.output',width = "50%"),br(), br(),br(), br(),downloadButton(outputId = "expdownload.output", label = "Download the expression of gene per sample"), br(),br(),br())),
            
            tabPanel(h4("Example files"),column(12, br(), br(), h4("Sample phenotype file: ",downloadButton("downloadmanifest", "Download sample phenotype file")), br(), br(), h4("Sample phenotype sheet with multi factors: ", downloadButton("downloadmanifestwithmfa", "Download sample phenotype file with multifactor columns")), br(), br(), h4("Raw counts file with gene IDs:",downloadButton("downloadrawcountsgeneid", "Download raw counts file with gene IDs")),br(), br(), h4("Raw counts file with gene symbols:",downloadButton("downloadrawcountsgenesymbol", "Download raw counts file with gene symbol"))))
                       # plotOutput('plot_sum3.output', width = "40%"),downloadButton(outputId = "sum2download.output", label = "Download the distribution of counts per sample")
            )))
            
            

