## iSeqQC: An Expression based Quality Control App  

### Introduction to iSeqQC
Quality Control (QC) in any high-throughput technology is a critical step, overlooking which can compromise the true conclusion of the data. As mentioned previously by others, a typical RNA-seq analysis should undergo two stages (four steps) of QC: 
1) Stage 1: Laboratory based- Step 1: RNA quality
2) Stage 2: Computational based- Step 2: Sequencing quality <br/>
                                 Step 3: Mapping quality <br/>
                                 Step 4: Batch effects 

Although, several attempts have been made to help researchers to circumvent the errors and biases in their sequencing experiment with the help of several tools which detects sequencing and mapping quality (FastQC, Picard, RseQC, RNA-SeQC, QC3 and so on). The scrutinizing of batch effects in any high-throughput sequencing experiment are often ignored or conducted on limited basis.  

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
  
### Input files requirement
iSeqQC requires two files for the analysis. Both files should be ASCII formatted **tab-delimited** file only
- File 1- Sample phenotype data: **Strictly limited to 4 columns in the same order as mentioned below (names case-sensitive)** </br>
**Sample names in first column 'samples' should match the names in counts matrix file**

column 1: samples<br/>
column 2: shortnames<br/>
column 3: groups<br/>
column 4: include<br/>

Example:<br/>

<table >
  <tr>
    <th>samples</th>
    <th>shortnames</th>
    <th>groups</th>  
    <th>include</th>  
  </tr>
  <tr>
    <td>Control_1</td>
    <td>C_1</td>  
    <td>control</td>  
    <td>TRUE</td>
  </tr>
  <tr>
    <td>Control_2</td>
    <td>C_2</td>  
    <td>control</td>  
    <td>TRUE</td>
  </tr>
  <tr>
    <td>Control_3</td>
    <td>C_3</td>  
    <td>control</td>  
    <td>TRUE</td>
  </tr>
  <tr>
    <td>Treated_1</td>
    <td>T_1</td>  
    <td>treated</td>  
    <td>TRUE</td>
  </tr>
  <tr>
    <td>Treated_2</td>
    <td>T_2</td>  
    <td>treated</td>  
    <td>TRUE</td>
  </tr>
  <tr>
    <td>Treated_3</td>
    <td>T_3</td>  
    <td>treated</td>  
    <td>TRUE</td>
  </tr>
</table>

- File 2- counts matrix file: **First column of this file should have official gene symbols under the name "gene_symbol"(case-sensitive)** 

Example:<br/>

<table >
  <tr>
    <th>gene_symbol</th>
    <th>Control_1</th>
    <th>Control_2</th>  
    <th>Control_3</th>
    <th>Treated_1</th>
    <th>Treated_2</th>
    <th>Treated_3</th>
  </tr>
 <tr>
    <td>TSPAN6</td>
    <td>642</td>
    <td>329</td>  
    <td>704</td>
    <td>507</td>
    <td>524</td>
    <td>629</td>
  </tr>
   <tr>
    <td>DPM1</td>
    <td>1443</td>
    <td>734</td>  
    <td>1502</td>
    <td>1175</td>
    <td>1543</td>
    <td>1111</td>
  </tr>
</table>

### Workflow
<img src= "https://user-images.githubusercontent.com/10853956/44124859-5ad6b14a-9ffd-11e8-9118-73f80c3b3514.png" width="500" height="250">


### Results Output
iSeqQC displays the results in a form of bunch of plots segregated into three different tabs in web UI: 
- Summary plots:
    - Summary Statistics of samples:  Provides summary statistics of the samples
    - Distribution of counts per sample: Plot displaying the distribution per sample (log2 transformed)
    - Density transcripts per sample: Plot displaying the density of distribution per sample (log2 transformed)
- Clustering plots:
  - Multi-Dimensional Scaling plot per sample: Plot displaying distances between each sample
  - Principal Component variances of all samples: Plot showing variations between each sample and segregation of groups
  - CPM variances per sample for the housekeeping genes: Plot showing variations in expression of housekeeping genes (log2 transformed)
- Correlational plots: 
  -  

