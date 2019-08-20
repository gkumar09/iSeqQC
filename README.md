## iSeqQC: An Expression based Quality Control tool

### Webpage
iSeqQC is readily available at:<br/>
<br/>
&nbsp;&nbsp;&nbsp;https://<span></span>http://cancerwebpa.jefferson.edu/iSeqQC/

### Command line
You can download iSeqQC from github using following (if git is already installed in your computer):<br/>
<br/>
&nbsp;&nbsp;&nbsp;git clone https://<span></span>github.com/gkumar09/iSeqQC.git

  
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


### Results Output
iSeqQC displays the results in a form of a summary table and several plots: 
    - Summary Statistics
    - Counts distribution
    - Mapped read density
    - Housekeeping gene expression
    - Principal Component variances
    - Hierarchical relationship between samples
    - Correlations (Pearson and Spearman)
  
### Workflow
<img src= "https://user-images.githubusercontent.com/10853956/44124859-5ad6b14a-9ffd-11e8-9118-73f80c3b3514.png" width="500" height="250">


