## iSeqQC: An Expression based Quality Control tool

It is an expression-based quality control tool to detect outliers either produced by batch effects or merely due to dissimilarity within a phenotypic group.
It can be utilized by three ways:

### Webpage (No Installation required)
iSeqQC is readily available at:<br/>
<br/>
&nbsp;&nbsp;&nbsp;<span></span>http://cancerwebpa.jefferson.edu/iSeqQC/

### Command line (Local R and Bioconductor libraries installation required)
Running iSeqQC locally requires: 
- Local installation of R or RStudio (version 3.5 or later)- if not available use https://cran.r-project.org/ to download.
- Installation of bioconductor packages using following commands: <br/>
	```
     if (!requireNamespace("BiocManager", quietly = TRUE))
     install.packages("BiocManager")
	 BiocManager::install()
     BiocManager::install(c("shiny", "FactoMineR", "factoextra", "som", "psych", "data.table", "ape", "corrplot", "limma", "DESeq2"))
    ```
- After successfully installing R/RStudio and related packages, iSeqQC can simply be run from 'iSeqQC_cli' directory using following command:<br/>

```
Rscript --vanilla iSeqQC_cli/iSeqQC.R exampleData/samplemanifestfile.txt {sample_phenotype_file} exampleData/genesymbol_rawcounts.txt {count_matrix} R {type_of_reads} SYMBOL {type_of_gene_identifier} H {Organism}
```

where,<br/>
type_of_reads: R for raw reads and N for normalized reads<br/>
type_of_gene_identifier: SYMBOL if count matrix has gene_symbols in first column and ID if it has gene_ids<br/>
Organism: H for Human, M for Mouse and O for others<br/>

### Local shiny installation (Local R and Bioconductor libraries installation required)
#### Prerequisities
Running iSeqQC locally requires: 
- Local installation of R or RStudio (version 3.5 or later)- if not available use https://cran.r-project.org/ to download.
- Installation of bioconductor packages using following commands: <br/>
	```
     if (!requireNamespace("BiocManager", quietly = TRUE))
     install.packages("BiocManager")
	 BiocManager::install()
     BiocManager::install(c("shiny", "FactoMineR", "factoextra", "som", "psych", "data.table", "ape", "corrplot", "limma", "DESeq2"))
    ```
- After successfully installing R/RStudio and related packages, iSeqQC can simply be run using following commands in R console:
```
setwd("path_to_local_iSeqQC_installation_directory")
library("shiny")
runApp("iSeqQC")
```

Please note: iSeqQC shiny installations are successful tested on safari v-12.1, chrome  v-79.0, and firefox v-72.2. However, we recommend google chrome for optimum usage. 


### Input files requirement
iSeqQC requires two files for the analysis. Both files should be ASCII formatted **tab-delimited** file only
- File 1- Sample phenotype data: **First 4 columns should strictly match the names and order as mentioned below (names case-sensitive)** </br>
**Sample names in first column 'samples' should match the names in counts matrix file**

column 1: samples<br/>
column 2: shortnames<br/>
column 3: groups<br/>
column 4: include<br/>
column 5-11: any factors such as library method, protocol etc.<br/>

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

- File 2- counts matrix file: **First column of this file should have official gene symbols or gene ids under the name "gene"(case-sensitive)** 

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

or

<table >
  <tr>
    <th>gene_id</th>
    <th>Control_1</th>
    <th>Control_2</th>  
    <th>Control_3</th>
    <th>Treated_1</th>
    <th>Treated_2</th>
    <th>Treated_3</th>
  </tr>
 <tr>
    <td>ENSG00000000003</td>
    <td>642</td>
    <td>329</td>  
    <td>704</td>
    <td>507</td>
    <td>524</td>
    <td>629</td>
  </tr>
   <tr>
    <td>ENSG00000000005</td>
    <td>1443</td>
    <td>734</td>  
    <td>1502</td>
    <td>1175</td>
    <td>1543</td>
    <td>1111</td>
  </tr>
</table>

### Results Output
iSeqQC displays the results in a form of a summary table and several plots: Summary statistics, counts distribution, Mapped read density, Housekeeping gene expression, Principal Component variances (zscored normalized), Principal Component variances (un-normalized), Hierarchical relationship between samples, Pearson correlation, Spearman correlation, GC bias)
  
### Workflow
<img src= "https://user-images.githubusercontent.com/10853956/44124859-5ad6b14a-9ffd-11e8-9118-73f80c3b3514.png" width="500" height="250">


