library("shiny")


shinyUI(fluidPage
        (theme = "bootstrap.css",
        titlePanel(h2("iSeqQC- An Expression based Quality Control tool", style= "font-family: 'American Typewriter'; color:#081d58"), windowTitle = "iSeqQC"),
        hr(),
          navbarPage("",
                     tabPanel(h4("About"), 
                              titlePanel(h3("Welcome to iSeqQC",align="center",style = "font-family: 'times'; color:#081d58; font-si16pt; line-height:1.8", br(), 
                                            h4("It is an expression-based quality control tool to detect outliers either produced by batch effects or merely due to dissimilarity within a phenotypic group.", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.8"),br(),
                                            h4("Input: It requires sample phenotype file and either raw counts matrix or normalized counts matrix (TPM/FPKM/RPKM etc.) file in ASCII formatted tab-delimited file as input. Type of file-type (raw or normalized counts), organism of study (human or mouse) and gene-identifier (gene id or symbol) is also required for successful execution of algorithm. Please make sure to match the column names in count matrix with 'samples' column in sample phenotype sheet. Also, please make sure that first 4 column names match the case and order as shown in the example files", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.8"),br(),
                                            h4("Output: It is presented in two tabs, where in ‘Quality Control’, you are required to input the above-mentioned information and obtain output is in form of an expression summary table and several QC plots such as density plots, PCA, hierchichal clustering, correlation plots. In ‘Expression Plots’ tab, you can choose your gene of interest to plot expression profile among different phenotypic groups. ", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.8"),br(),
                                            h4("Additional information can be obtained at https://github.com/gkumar09/iSeqQC", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.8"),br(),
                                            h4("Citation: Preprint", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.4"),
                                            h4("Kumar G, Ertel A, Feldman G, Kupper J, Fortina P (2019). iSeqQC: A Tool for Expression-Based Quality Control in RNA Sequencing. BioRxiv. https://doi.org/10.1101/768101", align="left",style = "font-family: 'times'; font-si14pt; line-height:1.2"),br(),
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
            
            

