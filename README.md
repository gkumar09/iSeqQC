## iSeqQC: An App to perform Quality Controls on expression based high-throughput data  

### Introduction to iSeqQC
Quality Control (QC) in any high-throughput technology is a critical step, overlooking which can compromise the true conclusion of the data. As mentioned previously by others, a typical RNA-seq analysis should undergo two stages (four steps) of QC: 
1) Stage 1: Laboratory based- Step 1: RNA quality
2) Stage 2: Computational based- Step 2: Sequencing quality <br/>
                                 Step 3: Mapping quality <br/>
                                 Step 4: Batch effects 

Although, several attempts have been made to help researchers to circumvent the errors and biases in their sequencing experiment by developing several tools which detects sequencing and mapping quality (FastQC, Picard, RseQC, RNA-SeQC, QC3 and so on). The scrutinizing of batch effects in any high-throughput sequencing experiment are often ignored or conducted on limited basis.  

iSeqQC is developed to analyze batch effects in a high-throughput study that could be due to laboratory conditions, reagent lots, personnel differences, different experiment times, biological samples and so on. It is a light-weight web UI tool that requires only sample phenotype data and count matrix file (from any counts reading tool e.g. HTseq, RSEM) to generate several plots indicating possible biasness in an experiment. It implements several statistical methods (unsupervised clustering, coorelation coffiecients, classical multidimensional scaling, mean-variance relationship). Not restricted to only RNA-seq, iSeqQC can be used with any of the high-throughput technology which results in counts matrix produced by overlap of reads with genes such as microarray, RNA-seq, ATAC-seq etc.    

### Download
You can download iSeqQC from github using following (if git is already installed in your computer):<br/>
<br/>
&nbsp;&nbsp;&nbsp;git clone https://<span></span>github.com/gkumar09/iSeqQC.git

### Prerequisities
Running iSeqQC locally requires: 
- Local installation of R or RStudio (version 3.5 or later)- if not available use https://cran.r-project.org/ to download.
- Installation of bioconductor packages using following commands: <br/>
     source("https://<span></span>bioconductor.org/biocLite.R") <br/>
     biocLite(c("shiny", "FactoMineR", "factoextra", "som", "psych", "data.table", "ape", "corrplot", "limma"))

### Running iSeqQC
After successfully installing R/RStudio and related packages, iSeqQC can simply be run using following commands in R console:<br/>
<br/>
setwd("path_to_local_iSeqQC_installation_directory") <br/>
library("shiny")<br/>
runApp("iSeqQC")
  
##

# Input files requirement
iSeqQC requires two files for the analysis. Both files should be ASCII formatted **tab-delimited** file only
1. File 1- Sample phenotype data: **Strictly limited to 4 columns in same order as mentioned below (names case-sensitive)** 
<br/>
column 1: samples<br/>
column 2: shortnames<br/>
column 3: groups<br/>
column 4: include
</br>
Example:
<br/>
<table border= "1px solid black">
  <tr>
    <th>samples</th>
    <th>shortnames</th>
    <th>groups</th>  
    <th>include</th>  
  </tr>
  <tr>
    <th>Control_1</th>
    <th>C_1</th>  
    <th>control</th>  
    <th>TRUE</th>
  </tr>
  <tr>
    <th>Control_2</th>
    <th>C_2</th>  
    <th>control</th>  
    <th>TRUE</th>
  </tr>
  <tr>
    <th>Control_3</th>
    <th>C_3</th>  
    <th>control</th>  
    <th>TRUE</th>
  </tr>
  <tr>
    <th>Treated_1</th>
    <th>T_1</th>  
    <th>treated</th>  
    <th>TRUE</th>
  </tr>
  <tr>
    <th>Treated_2</th>
    <th>T_2</th>  
    <th>treated</th>  
    <th>TRUE</th>
  </tr>
  <tr>
    <th>Treated_3</th>
    <th>T_3</th>  
    <th>treated</th>  
    <th>TRUE</th>
  </tr>
</table>
    




