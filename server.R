##############################
#                            #
# GSEA+: From the basics     # 
#                            #
##############################

##############################
# 0. Defining the functions
##############################

# Loading the required libraries
library(ggplot2)
library(shiny) 
library(shinydashboard)
library(knitr) 
library(data.table)
library(rintrojs) 
library(shinyWidgets)
library(rmarkdown) 
library(cowplot) 
library(dplyr)
library(devtools)


source("GSEAfunctions.R")



##########################################################
##########################################################
##########################################################

server <- shinyServer(function(input, output, session) {
    
    # Implementation of the IntroBox to guide the user (Dashboard/GSEA charts)
    steps <- reactive(
      data.frame(
        element = c("#panel1", "#help_search", "#summary", "#downloadTable", "#choices", "#downloadPlot", "#leading_edge"),
        intro = c(
          "The classical GSEA plot consists of 3 panels: (i) enrichment profile, (ii) the gene set members appearance, and (iii) the ranked list metric. 
                  All panels will be displayed when running GSEA with single inputs (i.e a single preranked list and gene set).
                  For multiple inputs just the first panel will be created. 
                  Customization changes are directily applied, there is no need to re-run the analysis.
                  This plot has a brush function associated to draw a zoomed version, when the option is activated",
          "This funcitonality allows to insert gene symbols that in the case to be present in the analyzed sets will be highlighted in the plot.",
          "This table will display the information of the associated to the brushed area in the plot.
                  The button below allow to obtain information of the values contained in each of the columns of the summary table.",
          "Retrieve a .csv file with the above table.",
          "The loaded files will be displayed here. 
                  This fucntionality allows to choose which of the inputs are used to run the GSEA. 
                  Several inputs can be selected.",
          "Retrieve a .pdf file with the image of the current displayed plot.",
          "Allow to extract the leading-edge subset in a .csv file. when running GSEA with single inputs."
        ))
    )
    
    observeEvent(input$help,
                 introjs(session, 
                         options = list(steps = steps(),
                                        "nextLabel" = "Next",
                                        "prevLabel" = "Previous",
                                        "skipLabel" = "Skip")
                 )
    )
    
    # To track the preranks and gene sets that are loaded by the user.
    values <- reactiveValues(inPre = NULL, inCand = NULL)
    
    # Loading the pre-ranked lists
    output$select_prerank_file <- renderUI({
      list(hr(),
           checkboxGroupInput("selection_prerank", label = "Loaded preranked lists", choices = input$gene_list$name))
    })
    
    # Loading the gene sets
    output$select_gene_set_file <- renderUI({
      list(hr(),
           checkboxGroupInput("selection_candidate", label = "Loaded gene sets", choices = input$candidate_list$name))
    })
    
    # Selecting the loaded pre-ranks
    output$choose_curve <- renderUI({
      list(hr(),
           checkboxGroupInput("selection_curve", label = "Choose the preranked to display", choices = input$selection_prerank))
    })
    
    # Selecting the loaded gene sets
    output$choose_set <- renderUI({
      list(hr(),
           checkboxGroupInput("selection_set", label = "Choose the gene set to display", choices = input$selection_candidate))
    })
    
    # Save the pathways of the input files to later read the data that they contain
    observeEvent(input$selection_prerank, {values$inPre <- input$gene_list$datapath})
    observeEvent(input$selection_candidate, {values$inCand <- input$candidate_list$datapath})
    
    # Displays technical information of the loaded files (Loading Data/Selection)
    output$genelistfiles <- renderTable(input$gene_list)
    output$candidatefiles <- renderTable(input$candidate_list)
    
    # Defines the input pre-ranked list -> MULTIPLE INPUTS
    gene_list_loaded <- eventReactive(input$selection_prerank, {
      
      multiple_gene_list <- c()
      
      for(file in values$inPre){
        ext <- tools::file_ext(file)
        req(file)
        # shinyFeedback::feedbackWarning(file, ext != "rnk", "Please upload a rnk file" )
        validate(need(ext == c("rnk", "csv", "tsv"), "Please upload a rnk or csv/tsv file"))
        
        if(ext == "rnk"){
          input_gene_list <- read.table(file)
        } else{
          input_gene_list <- read.csv(file)
        }
        
        # Defining the columns of the input file
        names <- input_gene_list$V1 # Save the 1st column of the ranked list into a variable
        values <- input_gene_list$V2 # Save the 2nd column of the ranked list into a variable
        
        data <- reading_inputs(names, values) # Generate data frame
        
        genes_data <- list(data)
        multiple_gene_list <- c(multiple_gene_list, genes_data)
      }
      
      multiple_gene_list
    })
    
    # Defines the input gene set -> MULTIPLE INPUTS
    candidate_list_loaded <- eventReactive(input$selection_candidate, {
      
      multiple_candidate <- c()
      
      # Add loaded gene sets
      if(is.null(values$inCand) == FALSE){
        for(file_cand in values$inCand){
          ext_cand <- tools::file_ext(file_cand)
          req(file_cand)
          validate(need(ext_cand == c("csv", "tsv", "txt"), "Please upload a csv/tsv or txt file"))
          
          if(ext_cand == "txt"){
            candidate <- list(readLines(file_cand))
          } else{
            if(input$header == 2){
              input_candidate_list <- read.csv(file_cand, sep = "", header = F)
            } else{
              input_candidate_list <- read.csv(file_cand, sep = "", header = T)
              colnames(input_candidate_list) <- "V1"
            }
            input_candidate_list <- unique(input_candidate_list$V1) # Select the gene symbols and Remove repeated elements
            candidate <- list(input_candidate_list)
          }
          
          multiple_candidate <- c(multiple_candidate, candidate)
        }
      }
      
      # Add human hallmark
      if(input$available_gene_set == 1){
        multiple_candidate <- c(multiple_candidate, hallmark_human)
      }
      
      # Add mouse hallmark
      if(input$available_gene_set == 2){
        multiple_candidate <- c(multiple_candidate, hallmark_mouse)
      }
      
      multiple_candidate
      
    })
    
    # Get the name of each of the selected sets (Dashboard/GSEA charts)
    candidate_list_names <- eventReactive(input$selection_candidate, {
      
      multiple_names <- c()
      
      if(is.null(values$inCand) == FALSE){
        multiple_names <- c(multiple_names, input$candidate_list$name)
      }
      
      if(input$available_gene_set == 1){
        multiple_names <- c(multiple_names, multiple_description)
      }
      
      if(input$available_gene_set == 2){
        multiple_names <- c(multiple_names, multiple_description)
      }
      
      multiple_names
      
    })
    
    # Get the name of each of the selected pre-ranks (Dashboard/GSEA charts)
    gene_list_names <- eventReactive(input$selection_prerank, {
      
      multiple_names <- c()
      
      if(is.null(values$inPre) == FALSE){
        multiple_names <-c(multiple_names, input$gene_list$name)
      }
      
      multiple_names
    })
    
    # To know the sets that have been selected to be displayed in the plot -> Returns a vector
    candidate_select <- eventReactive(input$display_options_gene_set, {
      
      draw <- c()
      
      for(i in input$display_options_gene_set){
        draw <- c(draw, which(candidate_list_names() == i))
      }
      
      selected_candidate <- c()
      
      for(i in draw){
        selected_candidate <- c(selected_candidate, candidate_list_loaded()[i])
      }
      
      selected_candidate
      
    })
    
    # To know the pre-ranks that have been selected to be displayed in the plot -> Returns a vector
    prerank_select <- eventReactive(input$display_options_prerank, {
      
      draw <- c()
      
      for(i in input$display_options_prerank){
        draw <- c(draw, which(gene_list_names() == i))
      }
      
      selected_prerank <- c()
      
      for(i in draw){
        selected_prerank <- c(selected_prerank, gene_list_loaded()[i])
      }
      
      selected_prerank
      
    })
    
    # Choose the curves to be displayed 
    output$choices <- renderUI({
      list(
        checkboxGroupInput("display_options_prerank", label = "Select prerank:", choices = input$selection_prerank, selected = input$selection_prerank[1]),
        checkboxGroupInput("display_options_gene_set", label = "Select gene sets:", choices = candidate_list_names(), selected = candidate_list_names()[1]))
    })
    
    # Obtain the gene symbols within the leading-edge subset
    leading_edge <- reactive({

      leading_edge_df <- c()
      max = 0 
      
      for(i in length(data_multi())){
        
        if(data_multi()[[i]]$MES_value > 0){ # Considering positive enrichment
          leading_edge_df <- c(leading_edge_df, list(data_multi()[[i]]$data[data_multi()[[i]]$data$ranking <= data_multi()$MES_rank]$names))
          if(length(data_multi()[[i]]$data[data_multi()[[i]]$data$ranking <= data_multi()$MES_rank]$names) > max){
            max = length(data_multi()[[i]]$data[data_multi()[[i]]$data$ranking <= data_multi()$MES_rank]$names)
          }
        } else{ # Considering negative enrichment
          leading_edge_df <- c(leading_edge_df, list(data_multi()[[i]]$data[data_multi()[[i]]$data$ranking >= data_multi()$MES_rank]$names))
          if(length(data_multi()[[i]]$data[data_multi()[[i]]$data$ranking >= data_multi()$MES_rank]$names) > max){
            max = length(data_multi()[[i]]$data[data_multi()[[i]]$data$ranking >= data_multi()$MES_rank]$names)
          }
        }
      }
      
      les <- data.frame(matrix(NA,    # Create empty data frame
                               nrow = max,
                               ncol = length(leading_edge_df)))
      
      for(i in 1:length(leading_edge_df)){
        les[i] <- c(leading_edge_df[[i]], rep(NA, max - length(leading_edge_df[[i]])))
      }
      
      les
      
    })
    
    # Contains the data of the GSEA execution
    data_multi <- reactive({
      data_multi <- c()
      count = 0
      
      for(data_frame in prerank_select()){ # Traverse all the input gene lists
        for(list in 1:length(candidate_select())){ # Traverse all the candidate gene lists
          count = count + 1
          updateProgressBar(session = session, id = "progress", value = count,
                            total = length(prerank_select()) * length(candidate_select()))
          
          current_candidate <- c(candidate_select()[list]) # Select a candidate gene list
          current_candidate <- current_candidate[[1[1]]] # Change the format
          
          GSEA <- OurGSEA_plus(data_frame, current_candidate) # Call to perform all GSEA functions
          
          GSEA_out_data <- list(GSEA)
          data_multi <- c(data_multi, GSEA_out_data)
        }
      }
      data_multi
    })
    
    # Collect the computed statistics
    important_statistics <- reactive({
      stat_values <- c(Comparison = NULL, Pvalue = NULL, NES = NULL, FDR_Qvalue = NULL,
                       TAGS = NULL, LIST = NULL, SIGNAL = NULL)
      
      count = 1 # Keep track of the curent comparison number
      
      for(data_frame in prerank_select()){
        for(list in 1:length(candidate_select())){
          
          current_candidate <- c(candidate_select()[list]) # Select a candidate gene list
          current_candidate <- current_candidate[[1[1]]] # Change the format
          
          GSEA <- OurGSEA_plus(data_frame, current_candidate) # Call to perform all GSEA functions
          
          nes_perm <- c()
          
          for(mes in permutation_multi()[[count]]$ES_null){
            if(mes > 0){
              NES <- mes / mean(permutation_multi()[[count]]$pos_null)
            } else{
              NES <- abs(mes) / mean(permutation_multi()[[count]]$neg_null)
            }
            nes_perm <- c(nes_perm, NES)
          }
          
          if(GSEA$MES_value > 0){ # Considering positive enrichment
            leading_edge <- GSEA$data[GSEA$data$ranking <= GSEA$MES_rank, ]
            observed_NES <- GSEA$MES_value / mean(permutation_multi()[[count]]$pos_null)
            nom_pval <- sum(permutation_multi()[[count]]$ES_null >= GSEA$MES_value) / 1000
            FDR <- (sum(nes_perm > observed_NES)) / 1000
          } else{ # Considering negative enrichment
            leading_edge <- GSEA$data[GSEA$data$ranking >= GSEA$MES_rank, ]
            observed_NES <- abs(GSEA$MES_value) / mean(permutation_multi()[[count]]$neg_null)
            nom_pval <- sum(permutation_multi()[[count]]$ES_null <= GSEA$MES_value) / 1000
            FDR <- (sum(nes_perm < observed_NES)) / 1000
          }
          
          t = (sum(with(leading_edge, indicator == 1)) / sum(with(GSEA$data, indicator == 1)))
          l = (nrow(leading_edge) / nrow(GSEA$data))
          s = t * (1 - l) * (nrow(GSEA$data) / (nrow(GSEA$data) - length(current_candidate)))
          
          stat_values$Comparison <- c(stat_values$Comparison, count)
          stat_values$Pvalue <- c(stat_values$Pvalue, nom_pval)
          stat_values$NES <- c(stat_values$NES, observed_NES)
          stat_values$FDR_Qvalue <- c(stat_values$FDR_Qvalue, FDR)
          stat_values$TAGS <- c(stat_values$TAGS, t * 100)
          stat_values$LIST <- c(stat_values$LIST, l * 100)
          stat_values$SIGNAL <- c(stat_values$SIGNAL, s * 100)
          
          count = count + 1
        }
      }
      
      as.data.frame(stat_values)
    })
    
    # Get table with the statistics (from null ES distribution) -> IT IS NOT DISPLAYED
    output$statistics <- renderTable({
      important_statistics()
    })
    
    # Draw null ES distribution
    permutation_multi <- reactive({
      
      null_dist <- c()
      
      count1 = 0
      permutation_num = 1000 # Defining the number of times we will make the random sampling
      total_count = length(prerank_select()) * length(candidate_select()) * 2
      
      for(i in 1:length(data_multi())){
        # Update the progress bar -> start of an iteration
        count1 = count1 + 1
        updateProgressBar(session = session, id = "progress_permutation", value = count1, 
                          total = total_count)
        
        ES_nulldist <- phenotype_permutation(data_multi()[[i]]$data, permutation_num) # Call to perform the phenotype permutation
        
        null_dist <- c(null_dist, list(ES_nulldist))
        
        # Update the progress bar -> end of an iteration 
        count1 = count1 + 1
        updateProgressBar(session = session, id = "progress_permutation", value = count1, 
                          total = total_count)
      }
      
      null_dist
    })
    
    # Collect the most important values 
    important_values <- reactive({
      summary_values <- c(Comparison = NULL, Prerank = NULL, GeneSet = NULL, NumPrerank = NULL,  NumGeneSet = NULL,
                          EnrichedClass = NULL, MES = NULL, MES_rank = NULL, NumLeading = NULL)
      
      count_prerank = 1
      count_gene_set = 1
      count_comparison = 1
      
      for(data_frame in prerank_select()){ # Traverse all the input gene lists
        for(list in 1:length(candidate_select())){ # Traverse all the candidate gene lists
          
          current_candidate <- c(candidate_select()[list]) # Select a candidate gene list
          current_candidate <- current_candidate[[1[1]]] # Change the format
          
          GSEA <- OurGSEA_plus(data_frame, current_candidate) # Call to perform all GSEA functions
          
          summary_values$Comparison <- c(summary_values$Comparison, count_comparison)
          summary_values$Prerank <- c(summary_values$Prerank, input$display_options_prerank[count_prerank])
          summary_values$GeneSet <- c(summary_values$GeneSet, input$display_options_gene_set[count_gene_set])
          summary_values$NumPrerank <- c(summary_values$NumPrerank, nrow(data_frame))
          summary_values$NumGeneSet <- c(summary_values$NumGeneSet, nrow(current_candidate))
          
          if(GSEA$MES_value > 0){
            summary_values$EnrichedClass <-c(summary_values$EnrichedClass, "Positive")
          } else{
            summary_values$EnrichedClass <- c(summary_values$EnrichedClass, "Negative")
          }
          
          summary_values$MES <- c(summary_values$MES, round(GSEA$MES_value, 5))
          summary_values$MES_rank <- c(summary_values$MES_rank, GSEA$MES_rank)
          # summary_values$NumOverlap <- c(summary_values$NumOverlap, length(which(GSEA$data$indicator == 1)))
          
          if(GSEA$MES_value > 0){
             summary_values$NumLeading <-c(summary_values$NumLeading, GSEA$MES_rank)
          } else{
             summary_values$NumLeading <- c(summary_values$NumLeading, nrow(GSEA$data) - GSEA$MES_rank)
          }
        
          count_gene_set = count_gene_set + 1
          count_comparison = count_comparison + 1
          
        }
        
        count_prerank = count_prerank + 1
        count_gene_set = 1 # Reset the count, start again for the new prerank
      }
      as.data.frame(summary_values)
    })
    
    # Get table with the most important values (Dashboard/ Significance)
    output$summary_copy <- renderTable({
      important_values()
    })
    
    # Get table with the most important values (Dashboard/ GSEA charts)
    output$summary <- renderTable({
      important_values()
    })
    
    # Get the classical GSEA plot -> MULTIPLE INPUTS (Displays the users' choices)
    draw_panel1 <- reactive({
      
      whole_data <- data.frame()
      count = 1
      
      for(data_frame in prerank_select()){ # Traverse all the input gene lists
        for(list in 1:length(candidate_select())){ # Traverse all the candidate gene lists
          
          current_candidate <- c(candidate_select()[list]) # Select a candidate gene list
          current_candidate <- current_candidate[[1[1]]] # Change the format
          
          GSEA <- OurGSEA_plus(data_frame, current_candidate) # Call to perform all GSEA functions
          
          group <- rep(count, nrow(GSEA$data)) # Create vector to show the category of each dataframe
          GSEA$data["comparison"] <- group
          
          # Update the dataframe, add each dataframe at the end each time
          whole_data <- rbind(whole_data, GSEA$data)
          
          count = count + 1
          
        }
      }
      
      GSEA_plot <- ggplot() +
        theme(
          panel.background = element_rect(fill = input$color_panel),
          plot.background = element_rect(fill = input$color_plot),
          plot.title = element_text(color = input$color_title, size = input$title_size, face = "bold", hjust = input$title_position)) + 
        geom_line(data = whole_data, mapping = aes(x = ranking, y = running_ES, color = factor(comparison)), size = input$curve_size) +
        labs(x = "Rank in Ordered Dataset", y = "Running Enrichment Score", title = input$description, color = "Comparison Number") 
      
      if(input$add_0 == 1){
        GSEA_plot <- GSEA_plot + geom_hline(mapping = aes(yintercept = 0), linetype = "dashed")
      }
      
      # When just one prerank and one gene set are selected, draw all the panels
      if(length(data_multi()) == 1){
        
        GSEA_plot <- GSEA_plot + theme(
          panel.background = element_rect(fill = input$color_panel),
          plot.background = element_rect(fill = input$color_plot),
          plot.title = element_text(color = input$color_title, size = input$title_size, face = "bold", hjust = input$title_position),
          legend.position = "none") + 
          geom_line(data = data_multi()[[1]]$data, mapping = aes(x = ranking, y = running_ES), color = input$color, size = input$curve_size) +
          labs(x = "Rank in Ordered Dataset", y = "Running Enrichment Score", title = input$description, color = "Comparison Number") 
        
        # GSEA_plot <- GSEA_plot + theme(legend.position = "none")
        
        # Add MES line
        if(input$add_mes == 1){
          GSEA_plot <- GSEA_plot + geom_vline(mapping = aes(xintercept = data_multi()[[1]]$MES_rank), linetype = "dashed", color = input$color_dashed)
        }
        
        # Add panel 2 (barcode + color range)
        if(input$add_panel2 == 1){
          up_colors <- colorRampPalette(c("red", "#FFFFFF")) # Ranges from red to white
          down_colors <- colorRampPalette(c("#FFFFFF", "#1b98e0")) # Ranges from white to blue
          
          hit_rank <- data.frame(which(data_multi()[[1]]$data$indicator == 1))
          indicator1 <- which(data_multi()[[1]]$data$indicator == 1)
          
          min <- min(data_multi()[[1]]$data$running_ES)
          max <- max(data_multi()[[1]]$data$running_ES)
          
          pos_val <- sum(data_multi()[[1]]$data$values > 0)
          if(pos_val > 0){
            limit <- which(data_multi()[[1]]$data$values == min(data_multi()[[1]]$data$values[data_multi()[[1]]$data$values > 0]))
          } else{
            limit <- 0
          }
          
          GSEA_plot <- GSEA_plot + geom_rug(data = hit_rank, mapping = aes(x = indicator1), alpha = 0.5,
                                            sides = "b", color = input$color_panel2, length = unit(0.1, "npc")) + ylim(min - 0.1, max) + 
            geom_rug(data = data_multi()[[1]]$data, mapping = aes(x = ranking),
                     color = c(up_colors(limit), down_colors(nrow(data_multi()[[1]]$data) - limit))) +
            labs(x = NULL)
          
          # Highlight line when searched gene is overlapping
          if(is.null(input$selection_genes) == FALSE){
            genes <- input$selection_genes
            filtered_data <- data_multi()[[1]]$data %>% filter(names %in% genes)
            # filtered_data <- filtered_data %>% filter(indicator == 1)
            miss_filtered_data <- filtered_data %>% filter(indicator == 0)
            hit_filtered_data <- filtered_data %>% filter(indicator == 1)
            
            # Draws a line to represent the position of the gene (if the gene is present, indicator = 1)
            for(miss in miss_filtered_data){
              GSEA_plot <- GSEA_plot + geom_point(data = miss_filtered_data, mapping = aes(x = ranking, y = running_ES)) + 
                geom_text(data = miss_filtered_data, mapping = aes(x = ranking, y = running_ES, label = names), size = 5, fontface = 2, vjust = -1) 
            }
            
            for(hit in hit_filtered_data){
              GSEA_plot <- GSEA_plot + geom_point(data = hit_filtered_data, mapping = aes(x = ranking, y = running_ES), color = input$color_search_gene) + 
                geom_text(data = hit_filtered_data, mapping = aes(x = ranking, y = running_ES, label = names), size = 5, fontface = 2, vjust = -1, color = input$color_search_gene) 
            }
            
            # Draws a line to represent the position of the gene (if the gene is present, indicator = 1)
            GSEA_plot <- GSEA_plot + geom_rug(data = hit_filtered_data, mapping = aes(x = ranking), sides = "b", length = unit(0.15, "npc"), color = input$color_search_gene, size = 2)
          }
        }
        
        # Display panel 3
        if(input$add_panel3 == 1){
          GSEA_plot <- GSEA_plot +
            theme(
              axis.ticks.x = element_blank(),
              axis.text.x = element_blank())
          
          threshold <- max(abs(min(data_multi()[[1]]$data$values)), abs(max(data_multi()[[1]]$data$values)))
          
          p2 <- ggplot(data = data_multi()[[1]]$data, mapping = aes(x = ranking, y = values)) + geom_col() + 
            labs(y = "Ranked List Metric", x = "Rank in Ordered Dataset") +
            theme_classic() + ylim(-threshold, threshold)
          
          if(input$line_exp == 1){
            p2 <- p2 + geom_vline(mapping = aes(xintercept = limit[1]), linetype = "dashed", color = input$color_exp_line) 
          }
          
          if(input$red_blue == 2){
            p2 <- p2 + geom_col(color = input$color_exp)
          } else{
            p2 <- p2 + geom_col(color = ifelse(data_multi()[[1]]$data$values > 0, 'red', 'blue'))
          }
          
          GSEA_plot <- plot_grid(GSEA_plot, p2, ncol = 1, nrow = 2, rel_heights = c(2,1))
          
        } else{
          GSEA_plot <- GSEA_plot + labs(x = "Rank in Ordered Dataset")
        }
        
        # Multiple inputs
        for(i in 1:length(data_multi())){
          # Search genes option
          if(is.null(input$selection_genes) == FALSE){
            genes <- input$selection_genes
            filtered_data <- data_multi()[[i]]$data %>% filter(names %in% genes)
            hit_filtered_data <- filtered_data %>% filter(indicator == 1)
            miss_filtered_data <- filtered_data %>% filter(indicator == 0)
            
            # Draws a line to represent the position of the gene (if the gene is present, indicator = 1)
            for(miss in miss_filtered_data){
              GSEA_plot <- GSEA_plot + geom_point(data = miss_filtered_data, mapping = aes(x = ranking, y = running_ES)) + 
                geom_text(data = miss_filtered_data, mapping = aes(x = ranking, y = running_ES, label = names), size = 5, fontface = 2, vjust = -1) 
            }
            
            for(hit in hit_filtered_data){
              GSEA_plot <- GSEA_plot + geom_point(data = hit_filtered_data, mapping = aes(x = ranking, y = running_ES), color = input$color_search_gene) + 
                geom_text(data = hit_filtered_data, mapping = aes(x = ranking, y = running_ES, label = names), size = 5, fontface = 2, vjust = -1, color = input$color_search_gene) 
            }
          }
          
          # Add MES line -> IT JUST DIAPLAY THE MES LINE OF THE LAST ITERATION
          # if(input$add_mes == 1){
          #   GSEA_plot <- GSEA_plot + geom_vline(mapping = aes(xintercept = data_multi()[[i]]$MES_rank), linetype = "dashed", color = i+1)
          # }
        }
        
      }
      
      GSEA_plot
      
    })
    
    # Get the classical GSEA plot with the statistical values
    draw_panel1_copy <- reactive({
      
      whole_data <- data.frame()
      count = 1
      NES = 0
      nom_pval = 0
      
      for(data_frame in prerank_select()){ # Traverse all the input gene lists
        for(list in 1:length(candidate_select())){ # Traverse all the candidate gene lists
          
          current_candidate <- c(candidate_select()[list]) # Select a candidate gene list
          current_candidate <- current_candidate[[1[1]]] # Change the format
          
          GSEA <- OurGSEA_plus(data_frame, current_candidate) # Call to perform all GSEA functions
          
          group <- rep(count, nrow(GSEA$data)) # Create vector to show the category of each dataframe
          GSEA$data["comparison"] <- group
          
          # Update the dataframe, add each dataframe at the end each time
          whole_data <- rbind(whole_data, GSEA$data)
          
          nes_perm <- c()
          
          for(mes in permutation_multi()[[count]]$ES_null){
            if(mes > 0){
              NES <- mes / mean(permutation_multi()[[count]]$pos_null)
            } else{
              NES <- abs(mes) / mean(permutation_multi()[[count]]$neg_null)
            }
            nes_perm <- c(nes_perm, NES)
          }
          
          if(GSEA$MES_value > 0){ # Considering positive enrichment
            leading_edge <- GSEA$data[GSEA$data$ranking <= GSEA$MES_rank, ]
            observed_NES <- GSEA$MES_value / mean(permutation_multi()[[count]]$pos_null)
            nom_pval <- sum(permutation_multi()[[count]]$ES_null >= GSEA$MES_value) / 1000
            FDR <- (sum(nes_perm > observed_NES)) / 1000
          } else{ # Considering negative enrichment
            leading_edge <- GSEA$data[GSEA$data$ranking >= GSEA$MES_rank, ]
            observed_NES <- abs(GSEA$MES_value) / mean(permutation_multi()[[count]]$neg_null)
            nom_pval <- sum(permutation_multi()[[count]]$ES_null <= GSEA$MES_value) / 1000
            FDR <- (sum(nes_perm < observed_NES)) / 1000
          }
          
          t = (sum(with(leading_edge, indicator == 1)) / sum(with(GSEA$data, indicator == 1)))
          l = (nrow(leading_edge) / nrow(GSEA$data))
          s = t * (1 - l) * (nrow(GSEA$data) / (nrow(GSEA$data) - length(current_candidate)))

          count = count + 1
          
        }
      }
      
      GSEA_plot <- ggplot() +
        theme(
          panel.background = element_rect(fill = input$color_panel),
          plot.background = element_rect(fill = input$color_plot),
          plot.title = element_text(color = input$color_title, size = input$title_size, face = "bold", hjust = input$title_position)) + 
        geom_line(data = whole_data, mapping = aes(x = ranking, y = running_ES, color = factor(comparison)), size = input$curve_size) +
        labs(x = "Rank in Ordered Dataset", y = "Running Enrichment Score", title = input$description, color = "Comparison Number") 
      
      # Add y-axis = 0
      if(input$add_0 == 1){
        GSEA_plot <- GSEA_plot + geom_hline(mapping = aes(yintercept = 0), linetype = "dashed")
      }
      
      # When just one prerank and one gene set are selected, draw all the panels
      if(length(data_multi()) == 1){
        
        GSEA_plot <- GSEA_plot + theme(
          panel.background = element_rect(fill = input$color_panel),
          plot.background = element_rect(fill = input$color_plot),
          plot.title = element_text(color = input$color_title, size = input$title_size, face = "bold", hjust = input$title_position),
          legend.position = "none") + 
          geom_line(data = data_multi()[[1]]$data, mapping = aes(x = ranking, y = running_ES), color = input$color, size = input$curve_size) +
          labs(x = "Rank in Ordered Dataset", y = "Running Enrichment Score", title = input$description, color = "Comparison Number") 
        
        # GSEA_plot <- GSEA_plot + theme(legend.position = "none")
        
        x_text = nrow(data_multi()[[1]]$data)
        y_text = round(max(data_multi()[[1]]$data$running_ES), 2)
        text1 = paste("NES: ", round(observed_NES, 2))
        text2 = paste("P-value: ", nom_pval)
        text3 = paste("FDR (q-value): ", FDR)
        
        GSEA_plot <- GSEA_plot + annotate(geom = "text", x = (x_text - 5000), y = y_text, label = text1, size = 5) +
          annotate(geom = "text", x = (x_text - 5000), y = (y_text - 0.05), label = text2, size = 5) +
          annotate(geom = "text", x = (x_text - 5000), y = (y_text - 0.10), label = text3, size = 5)
        
        
        min <- min(data_multi()[[1]]$data$running_ES)
        max <- max(data_multi()[[1]]$data$running_ES)
        
        # Add MES line
        if(input$add_mes == 1){
          GSEA_plot <- GSEA_plot + geom_vline(mapping = aes(xintercept = data_multi()[[1]]$MES_rank), linetype = "dashed", color = input$color_dashed)
        }
        
        # Add panel 2 (barcode + color range)
        if(input$add_panel2 == 1){
          up_colors <- colorRampPalette(c("red", "#FFFFFF")) # Ranges from red to white
          down_colors <- colorRampPalette(c("#FFFFFF", "#1b98e0")) # Ranges from white to blue
          
          hit_rank <- data.frame(which(data_multi()[[1]]$data$indicator == 1))
          indicator1 <- which(data_multi()[[1]]$data$indicator == 1)
          
          pos_val <- sum(data_multi()[[1]]$data$values > 0)
          if(pos_val > 0){
            limit <- which(data_multi()[[1]]$data$values == min(data_multi()[[1]]$data$values[data_multi()[[1]]$data$values > 0]))
          } else{
            limit <- 0
          }
          
          GSEA_plot <- GSEA_plot + geom_rug(data = hit_rank, mapping = aes(x = indicator1), 
                                            sides = "b", color = input$color_panel2, length = unit(0.1, "npc")) + ylim(min - 0.1, max) + 
            geom_rug(data = data_multi()[[1]]$data, mapping = aes(x = ranking), color = c(up_colors(limit), down_colors(nrow(data_multi()[[1]]$data) - limit))) +
            labs(x = NULL)
          
          # Highlight line when searched gene is overlapping
          if(is.null(input$selection_genes) == FALSE){
            genes <- input$selection_genes
            filtered_data <- data_multi()[[1]]$data %>% filter(names %in% genes)
            # filtered_data <- filtered_data %>% filter(indicator == 1)
            miss_filtered_data <- filtered_data %>% filter(indicator == 0)
            hit_filtered_data <- filtered_data %>% filter(indicator == 1)
            
            # Draws a line to represent the position of the gene (if the gene is present, indicator = 1)
            for(miss in miss_filtered_data){
              GSEA_plot <- GSEA_plot + geom_point(data = miss_filtered_data, mapping = aes(x = ranking, y = running_ES)) + 
                geom_text(data = miss_filtered_data, mapping = aes(x = ranking, y = running_ES, label = names), size = 5, fontface = 2, vjust = -1) 
            }
            
            for(hit in hit_filtered_data){
              GSEA_plot <- GSEA_plot + geom_point(data = hit_filtered_data, mapping = aes(x = ranking, y = running_ES), color = input$color_search_gene) + 
                geom_text(data = hit_filtered_data, mapping = aes(x = ranking, y = running_ES, label = names), size = 5, fontface = 2, vjust = -1, color = input$color_search_gene) 
            }
            
            # Draws a line to represent the position of the gene (if the gene is present, indicator = 1)
            GSEA_plot <- GSEA_plot + geom_rug(data = hit_filtered_data, mapping = aes(x = ranking), sides = "b", length = unit(0.15, "npc"), color = input$color_search_gene, size = 2)
          }
        }
        
        # Display panel 3
        if(input$add_panel3 == 1){
          GSEA_plot <- GSEA_plot +
            theme(
              axis.ticks.x = element_blank(),
              axis.text.x = element_blank())
          
          threshold <- max(abs(min(data_multi()[[1]]$data$values)), abs(max(data_multi()[[1]]$data$values)))
          
          p2 <- ggplot(data = data_multi()[[1]]$data, mapping = aes(x = ranking, y = values)) + geom_col() + 
            labs(y = "Ranked List Metric", x = "Rank in Ordered Dataset") +
            theme_classic() + ylim(-threshold, threshold)
          
          if(input$line_exp == 1){
            p2 <- p2 + geom_vline(mapping = aes(xintercept = limit[1]), linetype = "dashed", color = input$color_exp_line) 
          }
          
          if(input$red_blue == 2){
            p2 <- p2 + geom_col(color = input$color_exp)
          } else{
            p2 <- p2 + geom_col(color = ifelse(data_multi()[[1]]$data$values > 0, 'red', 'blue'))
          }
          
          GSEA_plot <- plot_grid(GSEA_plot, p2, ncol = 1, nrow = 2, rel_heights = c(2,1))
          
        } else{
          GSEA_plot <- GSEA_plot + labs(x = "Rank in Ordered Dataset")
        }
        
      }
      
      # Multiple comparisons
      for(i in 1:length(data_multi())){
        
        # x_text = nrow(data_multi()[[1]]$data)
        # y_text = round(max(abs(data_multi()[[1]]$data$running_ES)), 2)
        # text1 = paste("NES: ", round(NES, 2))
        # text2 = paste("P-value: ", nom_pval)
        # text3 = paste("FDR (q-value): ", FDR)
        
        # GSEA_plot <- GSEA_plot + annotate(geom = "text", x = (x_text - 5000), y = y_text, label = text1, size = 3, color = i+1) +
        # annotate(geom = "text", x = (x_text - 5000), y = (y_text - 0.05), label = text2, size = 3, color = i+1) +
        # annotate(geom = "text", x = (x_text - 5000), y = (y_text - 0.10), label = text3, size = 3, color = i+1)
        
        # Search genes option
        if(is.null(input$selection_genes) == FALSE){
          genes <- input$selection_genes
          filtered_data <- data_multi()[[i]]$data %>% filter(names %in% genes)
          hit_filtered_data <- filtered_data %>% filter(indicator == 1)
          miss_filtered_data <- filtered_data %>% filter(indicator == 0)
          
          # Draws a line to represent the position of the gene (if the gene is present, indicator = 1)
          for(miss in miss_filtered_data){
            GSEA_plot <- GSEA_plot + geom_point(data = miss_filtered_data, mapping = aes(x = ranking, y = running_ES)) + 
              geom_text(data = miss_filtered_data, mapping = aes(x = ranking, y = running_ES, label = names), size = 5, fontface = 2, vjust = -1) 
          }
          
          for(hit in hit_filtered_data){
            GSEA_plot <- GSEA_plot + geom_point(data = hit_filtered_data, mapping = aes(x = ranking, y = running_ES), color = input$color_search_gene) + 
              geom_text(data = hit_filtered_data, mapping = aes(x = ranking, y = running_ES, label = names), size = 5, fontface = 2, vjust = -1, color = input$color_search_gene) 
          }
        }
        
        # Add MES line -> IT JUST DIAPLAY THE MES LINE OF THE LAST ITERATION
        # if(input$add_mes == 1){
        #   GSEA_plot <- GSEA_plot + geom_vline(mapping = aes(xintercept = data_multi()[[i]]$MES_rank), linetype = "dashed", color = i+1)
        # }
      }
      
      GSEA_plot
      
    })
    
    # Draws the classical GSEA plot (Dashboard/GSEA charts)
    output$panel1 <- renderPlot({
      draw_panel1()
    })
    
    # Draws the classical GSEA plot (Dashboard/ Significance)
    output$panel1_copy <- renderPlot({
      draw_panel1_copy()
    })
    
    # Drawing the null distributions 
    output$null_dist <- renderUI({
      
      lapply(1:length(data_multi()), function(i){
        
        # Generate an identification for each of the plots
        id <- paste0("comparison_", i)
        plotOutput(outputId = id)
        
        output[[id]] <- renderPlot({
          
          # We do not use ggplot as it is more time consuming and does not have any major improvement
          
          # histogram <- hist(permutation_multi()[[i]])
          # plot(histogram, col = ifelse(histogram$mids > 0, 'red', 'blue'), xlab = "ES", ylab = "Number of random values", main = "Null distribution") 
          # abline(v = data_multi()[[i]]$MES_value, col = "orange", lwd = 8) # The vertical line is not displayed
          
          # Display the plot using ggplot()
          
          permutation <- permutation_multi()[[i]]$ES_null
          null <- data.frame(permutation)
          
          NULL_plot <- ggplot(data = null, mapping = aes(x = permutation)) + geom_histogram(bins = 30) + theme_classic() +
            geom_vline(data = data_multi()[i]$data, mapping = aes(xintercept = data_multi()[[i]]$MES_value), linetype = "dashed", color = input$color_dashed1) +
            labs(x = "ES", y = "P(ES)", title = "Random ES Distribution", caption = id) +
            theme(
              panel.background = element_rect(fill = input$color_panel1),
              plot.background = element_rect(fill = input$color_plot1),
              plot.title = element_text(color = input$color_title1, face = "bold", hjust = 0.5))
          
          NULL_plot
          
        })
        
      })
    })
    
    # Initialize the coordinates for the zoomed version of the plot
    ranges <- reactiveValues(x = NULL, y =  NULL)
    
    # Get the plot with the zoom selection
    output$zoom <- renderPlot({
      GSEA_plot <- draw_panel1() + coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)
      
      GSEA_plot
    })
    
    # Manages the selection of PANEL 1 to make the zoom
    observe({
      brush <- input$brush_action
      if(!is.null(brush)){
        ranges$x <- c(brush$xmin, brush$xmax)
        ranges$y <- c(brush$ymin, brush$ymax)
      } else{
        ranges$x <- NULL
        ranges$y <- NULL
      }
    })
    
    # Adds some comments on the genes of interest functionality
    output$help_search <- renderUI({
      stra <- paste("The gene name has to be entered strictly, be careful with capital letters.")
      strb <- paste("A highlighted vertical line will be displayed just if the gene is present in the gene set (single input).")
      strbb <- paste("The gene will be labelled in a different color just if the gene is present in the gene set (multiple input).")
      strc <- paste("Multiple inputs are allowed.")
      
      HTML(paste(stra, strb, strbb, strc, sep = '<br/>'))
    })
    
    # Displays info to interpret the data on the summary table
    output$help_table <- renderUI({
      a <- paste(strong("Col. 1: "), "Displays the label used to identificate each of the comparisons between a pre-ranked list and a gene set. It is linked to the legend of the GSEA plot.")
      b <- paste(strong("Col. 2|3: "), "Contains the name of the pre-ranked lists and the gene sets.")
      c <- paste(strong("Col. 4|5: "), "Total sum of genes in each of the pre-ranked lists and gene sets.")
      d <- paste(strong("Col. 6: "), "Indicates the up-regulated class. Positive if enrichment in the left, and negative otherwise.")
      e <- paste(strong("Col. 7: "), "Gives the values of the Maximum Enrichment Score (MES). The MES value indicates the maximum deviance of the ES from 0.")
      f <- paste(strong("Col. 8: "), "Displays the rank were the MES is located.")
      g <- paste(strong("Col. 10: "), "Total sum of genes that form the leading-edge subset. The leading-edge subset is defined by the MES.") 
      h <- paste("If positive enrichement we will consider those genes before the MES, and if negative enrichment the genes after the MES will be considered.")
      
      HTML(paste(a, b, c, d, f, g, h, sep = '<br/>'))
    })
    
    # Generates the button to download the leading_edge -> Creates a table (.csv) with the rank order and the gene symbols
    # CREATES AN EMPTY FILE
    output$leading_edge <- downloadHandler(
      filename = function(){
        paste("leading_edge.csv", sep = "")
      },
      content = function(file){
        write.csv(x = leading_edge(), file = file)
      }
    )
    
    # Generates the button to download the summary table
    output$downloadTable <- downloadHandler(
      filename = function(){
        paste("SummaryTable.csv", sep = "")
      },
      content = function(file){
        write.csv(x = important_values(), file = file, sep = "", col.names = TRUE)
      }
    )
    
    # Generates button to download the plot (Dashboard/GSEA charts)
    output$downloadPlot <- downloadHandler(
      filename = function(){
        paste0("GSEAplot.pdf")
      },
      content = function(file){
        pdf(file)
        print(draw_panel1())
        dev.off()
      })
    
    # Generates button to download the copied plot (Dashboard/Significance)
    output$sig_plot <- downloadHandler(
      filename = function(){
        paste0("GSEAplot_stats.pdf")
      },
      content = function(file){
        pdf(file)
        print(draw_panel1_copy())
        dev.off()
      })
    
    # Generates button to download the copied plot (Dashboard/Significance)
    output$sig_table <- downloadHandler(
      filename = function(){
        paste0("Table_significance.csv")
      },
      content = function(file){
        pdf(file)
        write.csv(x = important_statistics(), file = file, sep = "", col.names = TRUE)
        dev.off()
      })
    
    # Generate button to download the null distribution (Dashboard/Significance)
    output$sig_null <- downloadHandler(
      filename = function(){
        paste0("GSEA_null_dist.pdf")
      },
      content = function(file){
        pdf(file)
        for(i in 1:length(data_multi())){
          
          permutation <- permutation_multi()[[i]]$ES_null
          null <- data.frame(permutation)
          
          NULL_plot <- ggplot(data = null, mapping = aes(x = permutation)) + geom_histogram(bins = 30) + theme_classic() +
            geom_vline(data = data_multi()[i]$data, mapping = aes(xintercept = data_multi()[[i]]$MES_value), linetype = "dashed", color = input$color_dashed1) +
            labs(x = "ES", y = "P(ES)", title = "Random ES Distribution", caption = id) +
            theme(
              panel.background = element_rect(fill = input$color_panel1),
              plot.background = element_rect(fill = input$color_plot1),
              plot.title = element_text(color = input$color_title1, face = "bold", hjust = 0.5))
          
          print(NULL_plot)
        }
        dev.off()
      })
    
    # important_statistics
    
    # Generates button to download the report
    output$downloadData <- downloadHandler(
      filename = "reportGSEAplus.pdf",
      content = function(file){
        # Copy the report file to a temporary directory before processing
        tempReport <- file.path(tempdir(), "gsea+report.Rmd")
        file.copy("gsea+report.Rmd", tempReport, overwrite = TRUE)
        
        # Set up the parameters to pass to the document
        params <- list(
          # Parameters for the plots
          add_title = input$add_title,
          description = input$description,
          title_size = input$title_size,
          title_position = input$title_position,
          color_panel = input$color_panel,
          color_plot = input$color_plot,
          color_title = input$color_title,
          color = input$color,
          color_dashed = input$color_dashed,
          curve_size = input$curve_size,
          add_mes = input$add_mes,
          add_0 = input$add_0,
          add_panel2 = input$add_panel2,
          color_panel2 = input$color_panel2,
          color_search_gene = input$color_search_gene,
          add_panel3 = input$add_panel3,
          color_exp_line = input$color_exp_line,
          color_exp = input$color_exp,
          red_blue = input$red_blue,
          line_exp = input$line_exp,
          color_dashed1 = input$color_dashed1,
          
          # Data
          gene_list = gene_list_loaded(),
          candidate_list = candidate_list_loaded(),
          num_prerank = length(gene_list_names()),
          num_geneset = length(candidate_list_names()),
          names_prerank = gene_list_names(),
          names_geneset = candidate_list_names(),
          selection_genes = input$selection_genes,
          
          # Tables
          important_values = important_values(),
          important_statistics = important_statistics())
        
        # Knit the document passing in the params list
        rmarkdown::render(tempReport, output_file = file,
                          params = params,
                          envir = new.env(parent = globalenv()))
      })
    
  }
  
)