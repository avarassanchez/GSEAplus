# GSEA+

## Table of Contents
1. What is GSEA+?
2. Setup

## 1. What is GSEA+ ? 
GSEA+ is an entirely new application that aim to revitalise the current Gene Set Enrichment Analysis (GSEA) resources, allowing interactive but also batch processing, to cope withthe challenges posed by community of GSEA users. The software is presented into two distinct modes: web page and command line.

Selection of most relevant GSEA+ functionalities:
* Multiplots of enrichments over several pre-ranked lists and gene sets of RNA-seq experiments.
* Huge range of customization options to modify which outputs are displayed and how they look alike.
* Obtain personalised ready-to-publish figures.
* Gene finder to locate (if possible) in the enrichment progfile the searched genes.
* Zoomed version of the classical GSEA plot to get a more precise view of the desired section.
* Generation of a compact and well-structured report to present the most relevant output and additional interpretations.

## 2. Setup

GSEA+ source code is entirely written in R. The latest full distribution release can be downloaded from GitHub:

https://github.com/avarassanchez/GSEAplus

List of available archives and folders:
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
