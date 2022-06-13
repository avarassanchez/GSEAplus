##############################
#                            #
# GSEA+: From the basics     # 
#                            #
##############################

##############################
# 0. Defining the functions
##############################

# Loading the required libraries
library(ggplot2) # Library to draw the plots
library(shiny) # Basic library to make the application work
library(dplyr)
library(shinydashboard) # Library to define a structure to the interface
library(devtools)
library(shinyFeedback)
library(progress)
library(rmarkdown) # Library to create the report
library(knitr) # Library to create the report
library(data.table)
library(cowplot) 
library(rintrojs) # Library to create the introductory help to the users
library(shinyWidgets)



source("GSEAfunctions.R")



##########################################################
##########################################################
##########################################################

ui <- shinyUI(dashboardPage(
  
  skin = "black",
  title = "GSEA+", # Application title
  
  # Application Header
  dashboardHeader(title = strong("GSEA+"), titleWidth = 230,
                  dropdownMenu(type = "notifications", headerText = strong("HELP"), icon = icon("question"), badgeStatus = NULL,
                               notificationItem(text = "Main panel where the results will be displayed", icon = icon("dashboard")),
                               notificationItem(text = "Allows to upload and select multiple gene lists (rnk) and candidate lists (csv/tsv).", icon = icon("spinner")),
                               notificationItem(text = "Enables the user to configure how graphs will be displayed.", icon = icon("chart-area")),
                               notificationItem(text = "Export the outcome of the analysis.", icon = icon("download")),
                               notificationItem(text = strong("All changes in the selection criteria need to be confirmed clicking in the confirm button.", icon = icon("alert", lib = "glyphicon")))),
                  tags$li(a(strong("About GSEA+"), height = 40, href = "Add link to GitHub", title = "", target = "_blank"),
                          class = "dropdown")),
  
  # Application sidebar
  dashboardSidebar(
    introjsUI(), # Initialize the IntroBox
    sidebarMenu(
      actionButton("run_gsea", label = strong("RUN GSEA"), icon = icon("play-circle"), width = 200),
      menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
      menuItem("Loading Data", tabName = "loading_data", icon = icon("spinner")),
      menuItem("Plots Customization", tabName = "customization", icon = icon("chart-area")),
      menuItem("Download Report", tabName = "download", icon = icon("download"),
               div(downloadButton(outputId = "downloadData", label = "GSEA+ Report", icon = icon("download"), style = "color: black; margin-left: 15px; margin-bottom: 5px;"))),
      HTML("<br><br><br>"),
      HTML("GSEA charts"),
      HTML("<br><br>"),
      progressBar(id = "progress", value = 0, display_pct = TRUE),
      HTML("Permutation method"),
      HTML("<br><br>"),
      progressBar(id = "progress_permutation", value = 0, display_pct = TRUE),
      HTML("<br><br><br>")
    )),
  
  # Application body
  dashboardBody(
    introjsUI(), # Initialize the IntroBox
    tabItems(
      tabItem(tabName = "dashboard",
              tabsetPanel(
                tabPanel(title = strong("PIPELINE"), icon = icon("random"),
                         fluidPage(
                           fluidRow(
                             column(width = 12, box(width = NULL, br(),
                                                    tags$strong("STEP 1: Assessing the overlap"),
                                                    tags$p("In this step we attempt to compare the gene and the candidate lists by defining two sets:
                                genes in both lists (hit) and genes not present in the candidate list (misses). To do so, we traverse all the genes
                                in the preranked gene set while updating a vector with those gene symbols that are shared between sets. We track the presence of 
                                the genes with a binomial variable: 1 if present and 0 if absent. This new variable will be added to the initial data frame.
                                Additionally count the number of genes that overlap between sets"),
                                                    tags$strong("STEP 2: Computing the Sk"),
                                                    tags$p("It is necessary to define a Sk value for those overlapping and non-overlapping genes. We will use the same approach as in the previous
                                 step, traversing the data frame and updating a new variable. The Sk value for overlapping genes depends on the correlation value (rank metric).
                                 The Sk value for non-overlapping genes is independent of the correlation value, thus, is a constant that just considers the number of unhits"),
                                                    tags$strong("STEP 3: Computing the ES"),
                                                    tags$p("To correctly compute the Enrichment Score (ES) we need to consider, as stated by Subramanian et al., the magnitude of the increment depends 
                                 on the correlation of the gene with the phenotype. The procedure is the same one: we create a vector to store the running sum and use a loop
                                 to traverse all the genes. Based on the presence variable (Step 1), we will add or subtract the Sk value (Step 2)."),
                                                    tags$strong("STEP 4: Drawing the GSEA plot"),
                                                    tags$p("The classical GSEA plot is formed by three panels. The top panel corresponds to the enriched profile, and the shape of the ES 
                                 relates to the arrangement of the genes in the gene set compared with the whole expression set. The higher peak matches the maximum ES. 
                                 The middle panel shows where the gene set members (i.e. genes annotated inside this ontology) appear in the ranked gene list. 
                                 The bottom panel represents the ranked list metric in a curve, which measures the gene’s correlation with a phenotype."),
                                                    tags$strong("STEP 5: Significance testing"),
                                                    tags$p("GSEA employs permutation methods (= resampling) to generate a null distribution for each gene set.
                                 We used a phenotype permutation to randomly swap the sample labels and recalculate each time the ES.
                                 With such approach we obtain how widely the ES varies and how often two groups are effectively the same.
                                 In our code we complete 1000 permutations in which we randomly assign the original phenotype labels, re-ordering the genes
                                 and re-computing the ES. Once completed, we create a histogram of the corresponding ES null distribution."),
                                                    tags$p("Compute the p-value"),
                                                    tags$p("Multiple testing correction (FDR)"))
                             )))),
                tabPanel(title = strong("GSEA CHARTS"), icon = icon("chart-area"),
                         column(width = 4, verticalLayout(box(plotOutput("panel1", brush = brushOpts(id = "brush_action", resetOnNew = TRUE)), # The output plot is related by the brush_action
                                                              div(materialSwitch(inputId = "on_off", label = "Display zoom", status = "danger", value = FALSE),
                                                                  downloadButton(outputId = "downloadPlot", label = "", icon = icon("download"), style = "color: black; margin-left: 15px; margin-bottom: 5px;"), 
                                                                  downloadButton(outputId = "leading_edge", label = "", icon = icon("barcode", lib = "glyphicon"), style = "color: black; margin-left: 15px; margin-bottom: 5px;"),
                                                                  actionButton("help", "", icon = icon("question"), style = "color: black; margin-left: 15px; margin-bottom: 5px;")),
                                                              width = 300),
                                                          conditionalPanel(condition = "input.on_off == 1",
                                                                           box(plotOutput("zoom"), width = 300)))), 
                         column(width = 8, verticalLayout(box(h4("Genes of interest"),
                                                              selectizeInput("selection_genes", label = "Introduce a gene symbol:",
                                                                             selected = NULL, multiple = TRUE, choices = c(),
                                                                             options = list("plugins" = list("remove_button"), "create" = TRUE, "persist" = TRUE)), 
                                                              htmlOutput("help_search"), width = 500),
                                                          box(h4("GSEA results summary"),
                                                              tableOutput("summary"), 
                                                              div(materialSwitch(inputId = "help_on_off", label = "Column guide", status = "danger", value = FALSE),
                                                                  downloadButton(outputId = "downloadTable", label = "", icon = icon("download"), style = "color: black; margin-left: 15px; margin-bottom: 5px;"),
                                                                  HTML("<br><br><br>"),
                                                                  conditionalPanel(condition = "input.help_on_off == 1",
                                                                                   htmlOutput("help_table"))),
                                                              width = 500),
                                                          box(h4("Select inputs to be displayed"),
                                                              p("At least a preranked list and a gene set needs to be selected."),
                                                              uiOutput("choices"), width = 500)
                         ))),
                tabPanel(title = strong("SIGNIFICANCE"), icon = icon("sync"),
                         column(width = 4, verticalLayout(box(plotOutput("panel1_copy"), width = 300),
                                                          box(uiOutput("null_dist"), width = 300))),
                         column(width = 8, verticalLayout(box(tableOutput("summary_copy"),
                                                              tableOutput("statistics"),
                                                              width = 500)))))),
      tabItem(tabName = "loading_data",
              tabsetPanel(tabPanel(strong("SELECTION"), icon = icon("list"),
                                   box(h2("Loading files:"),
                                       textOutput("count"),
                                       p("Remember that the preranked list needs to be in .rnk format and the gene set file in .tsv/.csv or .txt format."),
                                       fileInput("gene_list", label = "Load the preranked list", multiple = TRUE, accept = c(".rnk", ".csv", ".tsv"), placeholder = "No file selected"),
                                       tableOutput("genelistfiles"),
                                       fileInput("candidate_list", label = "Load the gene set", multiple = TRUE, accept = c(".tsv", ".csv", ".txt"), placeholder = "No file selected"),
                                       tableOutput("candidatefiles")),
                                   box(h2("Choosing files:"),
                                       p("At least a preranked list and a gene set must be selected to properly run GSEA. Multiple input selection is allowed."),
                                       uiOutput("select_prerank_file"),
                                       uiOutput("select_gene_set_file"),
                                       checkboxGroupInput("available_gene_set", label = "Available gene sets", choices = list("Hallmark Human (May 2022)" = 1, "Hallmark Mouse (May 2022)" = 2, "None" = 3), selected = 3)))
              )),
      tabItem(tabName = "customization",
              tabsetPanel(tabPanel(strong("GSEA CHARTS"), icon = icon("chart-area"),
                                   box(h2("Enrichement profile:"),
                                       h4("Title related modifications"),
                                       radioButtons("add_title", label = "Display title", choices = list("Yes" = 1, "No" = 2), selected = 2),
                                       conditionalPanel(condition = "input.add_title == 1",
                                                        textInput("description", label = "Specify a title", ""),
                                                        numericInput("title_size", label = "Letter's size", value = 15, min = 10, max = 100),
                                                        selectInput("title_position", label = "Position",
                                                                    choices = list("Left" = 0, "Center" = 0.5, "Right" = 1), selected = 0.5)),
                                       # sliderInput("ranks", label = "Range of ranks", min = 1, max = nrow(input$gene_list), value = c(0, nrow(input$gene_list))),
                                       h4("Modifying the colors"),
                                       selectInput("color_panel", label = "Panel background", 
                                                   choices = list("White" = 0, "Black" = 1, "Green" = 3, "Yellow" = 7, "Red" = 2, "Blue" = 4, "Cyan" = 5, "Magenta" = 6, "Gray" = 8), selected = 0),
                                       p("Refers to the space within the axis. Corresponds to the area where the curve is displayed."),
                                       selectInput("color_plot", label = "Plot background", 
                                                   choices = list("White" = 0, "Black" = 1, "Green" = 3, "Yellow" = 7, "Red" = 2, "Blue" = 4, "Cyan" = 5, "Magenta" = 6, "Gray" = 8), selected = 0),
                                       p("Refers to the space outside the axis. The name of the axis are displayed in this area."),
                                       conditionalPanel(condition = "input.add_title == 1", 
                                                        selectInput("color_title", label = "Title", 
                                                                    choices = list("Black" = 1, "Green" = 3, "Yellow" = 7, "Red" = 2, "Blue" = 4, "Cyan" = 5, "Magenta" = 6, "Gray" = 8), selected = 1)),
                                       selectInput("color", label = "Curve", 
                                                   choices = list("Black" = 1, "Green" = 3, "Yellow" = 7, "Red" = 2, "Blue" = 4, "Cyan" = 5, "Magenta" = 6, "Gray" = 8), selected = 3),
                                       conditionalPanel(condition = "input.add_mes == 1",
                                                        selectInput("color_dashed", label = "MES line", 
                                                                    choices = list("Black" = 1, "Green" = 3, "Yellow" = 7, "Red" = 2, "Blue" = 4, "Cyan" = 5, "Magenta" = 6, "Gray" = 8), selected = 2)),
                                       h4("Other parameters"),
                                       sliderInput("curve_size", label = "Curve's thickness", min = 0.5, max = 2, value = 1, step = 0.1),
                                       p("As the thickness increases, the curve's resolution gets reduced.")
                                   ),
                                   box(radioButtons("add_mes", label = "Add maximum ES", choices = list("Yes" = 1, "No" = 2), selected = 1),
                                       radioButtons("add_0", label = "Display ES = 0 line", choices = list("Yes" = 1, "No" = 2), selected = 1),
                                       h2("Gene set members appearance:"),
                                       radioButtons("add_panel2", label = "Display this panel", choices = list("Yes" = 1, "No" = 2), selected = 1),
                                       conditionalPanel(condition = "input.add_panel2 == 1",
                                                        selectInput("color_panel2", label = "Bar's color", 
                                                                    choices = list("Black" = 1, "Green" = 3, "Yellow" = 7, "Red" = 2, "Blue" = 4, "Cyan" = 5, "Magenta" = 6, "Gray" = 8), selected = 1),
                                                        selectInput("color_search_gene", label = "Highlighted line's color",
                                                                    choices = list("Orange" = "orange", "Green" = 3, "Yellow" = 7, "Red" = 2, "Blue" = 4, "Cyan" = 5, "Magenta" = 6, "Gray" = 8), selected = "orange")),
                                       h2("Ranked list metric:"),
                                       radioButtons("add_panel3", label = "Display this panel", choices = list("Yes" = 1, "No" = 2), selected = 1),
                                       conditionalPanel(condition = "input.add_panel3 == 1",
                                                        radioButtons("line_exp", label = "Add threshold of expression", choices = list("Yes" = 1, "No" = 2), selected = 2),
                                                        h4("Modifying colors"),
                                                        conditionalPanel(condition = "input.line_exp == 1",
                                                                         selectInput("color_exp_line", label = "Threshold line's color", 
                                                                                     choices = list("Black" = 1, "Green" = 3, "Yellow" = 7, "Red" = 2, "Blue" = 4, "Cyan" = 5, "Magenta" = 6, "Gray" = 8), selected = 6)),
                                                        radioButtons("red_blue", label = "Paint by expression", choices = list("Yes" = 1, "No" = 2), selected = 2),
                                                        conditionalPanel(condition = "input.red_blue == 2", 
                                                                         selectInput("color_exp", label = "Plain color", 
                                                                                     choices = list("Black" = 1, "Green" = 3, "Yellow" = 7, "Red" = 2, "Blue" = 4, "Cyan" = 5, "Magenta" = 6, "Gray" = 8), selected = 1))))
              ),
              tabPanel(strong("SIGNIFICANCE"), icon = icon("sync"),
                       box(h2("Customization options:"),
                           selectInput("color_exp1", label = "Plain color", 
                                       choices = list("Black" = 1, "Green" = 3, "Yellow" = 7, "Red" = 2, "Blue" = 4, "Cyan" = 5, "Magenta" = 6, "Gray" = 8), selected = 1),
                           selectInput("color_panel1", label = "Panel background", 
                                       choices = list("White" = 0, "Black" = 1, "Green" = 3, "Yellow" = 7, "Red" = 2, "Blue" = 4, "Cyan" = 5, "Magenta" = 6, "Gray" = 8), selected = 0),
                           p("Refers to the space within the axis. Corresponds to the area where the curve is displayed."),
                           selectInput("color_plot1", label = "Plot background", 
                                       choices = list("White" = 0, "Black" = 1, "Green" = 3, "Yellow" = 7, "Red" = 2, "Blue" = 4, "Cyan" = 5, "Magenta" = 6, "Gray" = 8), selected = 0),
                           p("Refers to the space outside the axis. The name of the axis are displayed in this area."),
                           conditionalPanel(condition = "input.add_title1 == 1", 
                                            selectInput("color_title1", label = "Title", 
                                                        choices = list("Black" = 1, "Green" = 3, "Yellow" = 7, "Red" = 2, "Blue" = 4, "Cyan" = 5, "Magenta" = 6, "Gray" = 8), selected = 1)),
                           selectInput("color_dashed1", label = "Observed MES line", 
                                       choices = list("Black" = 1, "Green" = 3, "Yellow" = 7, "Red" = 2, "Blue" = 4, "Cyan" = 5, "Magenta" = 6, "Gray" = 8), selected = 2),
                           p("The observed MES corresponds to the maximum enrichment score computed from the original data (without any permutation method)."))
              ))
      )
    )
  )
  
)) 