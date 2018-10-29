library("shiny")
#library("FactoMineR")

shinyUI(fluidPage
        (theme = "bootstrap.css",
        headerPanel(h2("iSeqQC- An Expression based Quality Control app", style= "font-family: 'American Typewriter'; color:#081d58")),
        hr(),
        sidebarPanel(
          fileInput("manifest", label="Choose sample phenotype File", accept= "text/plain"), helpText("Please remember to use correct column names and matched sample names to avoid any errors")),
        sidebarPanel(
          tags$head(tags$style(type="text/css", "#loadmessage {position: fixed;top: 350px;left: 0px;width: 100%;padding: 5px 0px 5px 0px;text-align: center;font-weight: bold;font-size: 100%;color: #000000; z-index: 105;}"),
                    conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                     tags$div("Please wait while the server is working...",id="loadmessage"))),
          tags$head(
            tags$style(HTML(" .shiny-output-error-validation {color: red;}"))
            ),
          fileInput("countsmatrix", label="Choose Counts matrix File", accept= "text/plain"), helpText("Please remember to have official gene symbols in the 'gene_symbol' column (not IDs) to avoid any errors")),
          tabsetPanel(type="tabs",
            tabPanel("Summary module", column(12, align="center", br(), verbatimTextOutput(outputId= 'table1.output'), br(), br(), plotOutput('plot_sum2.output', width = "40%"),downloadButton(outputId = "sum1download.output", label = "Download the distribution of counts per sample"), br(), br(), plotOutput('plot_sum3.output', width = "50%"),downloadButton(outputId = "sum2download.output", label = "Download the detected transcripts per sample"),br(), br(),br())),
            tabPanel("Clustering module", column(12, align="center", br(), br(), plotOutput('plot2.output', width = "80%"),downloadButton(outputId = "pca2download.output", label = "Download the Principal components plot for all samples"), br(), br(),plotOutput('plot4.output', width = "80%", height = "400px"),downloadButton(outputId = "hierdownload.output", label = "Download the hierarchichal plot"), br(), br(), plotOutput('plot5.output', width = "80%", height = "500px"),downloadButton(outputId = "cordownload.output", label = "Download the pearson correlation plot"), br(), br(),br()))
            #tabPanel("Correlational",column(12, align="center", br(), br(), plotOutput('plot4.output', width = "80%", height = "400px"),downloadButton(outputId = "hierdownload.output", label = "Download the hierarchichal plot"), br(), br(), plotOutput('plot5.output', width = "80%", height = "500px"),downloadButton(outputId = "cordownload.output", label = "Download the pearson correlation plot"), br(), br(),br() ))
            ))
)

