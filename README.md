# GSEA+

## Table of Contents
1. What is GSEA+ ?
2. Setup

## 1. What is GSEA+ ? 
GSEA+ is an entirely new application that aim to revitalise the current Gene Set Enrichment Analysis (GSEA) resources, allowing interactive but also batch processing, to cope with the challenges posed by community of GSEA users. We plan to improve the current GSEA application by introducing novel functions that complement the ones now available and confer higher flexibility to the current GSEA pipelines: (i) to produce an interactive interface; (ii) to enable multiple ranking comparisons on the same set or on several RNA-seq groups; (iii) to facilitate the interpretation of the results; and (iv) to create personalised high-resolution and ready-to-publish plots.

The software is presented into two distinct modes: dynammic application and integrated in command-line pipelines.

Selection of most relevant GSEA+ functionalities:
* Multiplots of enrichments over several pre-ranked lists and gene sets of RNA-seq experiments.
* Huge range of customization options to modify which outputs are displayed and how they look alike.
* Obtain personalised ready-to-publish figures.
* Gene finder to locate (if possible) in the enrichment progfile the searched genes.
* Zoomed version of the classical GSEA plot to get a more precise view of the desired section.
* Generation of a compact and well-structured report to present the most relevant output and additional interpretations.

## 2. Setup

GSEA+ source code is entirely written in R. The latest full distribution release can be downloaded from GitHub:
```
https://github.com/avarassanchez/GSEAplus
```
### List of Files and Folders
* README: Description of the project and guidelines to run GSEA+.
* server.R: Server function definition of R-Shiny application.
* ui.R: Interface definition of R-Shiny application.
* GSEAfunctions: Definition of the functions required to run a GSEA analysis. 
* app.R: Call to the two-files R-Shiny application.
* /Data:
  * ES_NP_preranked_paired.rnk: Pre-ranked file to use as part of the demo.
  * QuickGO_Ossification.csv: Gene set file to use as part of the demo.
  * h.all.v7.5.1.symbols.gmt: Hallmark gene set collection from MSigDB (May 2022).
  * mart_export.txt: List of orthologous human-mouse nomenclature translator.

### How to Run GSEA+ 
GSEA+ can be executed locally from R-Studio, by command-line or using a browser.

The required R libraries are listed below:
```
library(ggplot2)
library(shiny) 
library(shinydashboard)
library(knitr) 
library(data.table)
library(rintrojs) 
library(shinyWidgets)
library(rmarkdown) 
```
Executing GSEA+ using R-Studio is very straight format: open the app.R and click on Run App. To compile GSEA+ through the terminal it is necessary to access the GSEAplus directory and execute the command below. If no errors, a link will be displayed with the same format as: "Listening on http:// XXX.X.X.X:YYYY"
```
Rscript app.R
```

### List of Commands
| Functions | Description | Input | Output | 
| --------- | ----------- | ----- | ------ |
| reading_inputs() | Filters the pre-ranked input lists. The resulting data frame consists of the two columns from the input file and a ranking variable to index the genes | Columns of the pre-ranked list | Ordered data frame |
| assessing_overlapping() | Makes use of a binomial variable to assess if a gene is present in the pre-ranked list and in the gene set (overlaps) | Data frame from reading_inputs() function | Updated data frame |
| computing_sk_ES() | Computes the Enrichment Score (ES) | Data frame from assessing_overlap() and number of overlapping genes | Updated data frame |
| phenotype_permutation() | Use of random permutation methods to generate a null distribution to assess the statistical significance of the ES | Data frame from computing_sk_ES() and number of permutations to be completed | Vector with the stored Maximum Enrichment Scores (MES) of each permutation |
| OurGSEAplus() | Main function call to execute GSEA | Original data frame from reading_inputs() and the input gene set | The final updated data frame from computing_sk_ES() and the MES rank (order) and the MES value (score) |


### Inputs Data Formats
2 input files are necessary to execute GSEA: pre-ranked file and gene set. 
* The pre-ranked list are obtained from the RNA-seq pipeline. Consist of a two-column format containing the gene symbol and the gene rank metric (from differential expression analysis with softwares such as DESeq2). 
* The gene set consist of a gene symbol lists related to a particular category. Can be loaded from two different sources: users' own signatures or an already available gene set collection. These approaches enable a more precise or exploratory search, respectively. The Hallmarks collection from The Molecular Signatures Database (MSigDB) has been included containing a total of 50 gene sets. 
