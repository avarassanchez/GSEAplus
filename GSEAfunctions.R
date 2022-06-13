##########################################################
##########################################################
##########################################################

# Functions

# Takes the names and values of the preranked list (1st and 2nd column)
# Returns the dataframe with the rank order
reading_inputs <- function(names, values){
  
  # Creating a data frame
  data <- data.frame(
    names,
    values
  )
  
  data <- data[order(data$values),] # Ordering the data frame by the value of expression
  
  ranking <- c(length(data$names):1)
  data["ranking"] <- ranking
  
  data <- data[order(data$ranking),] # Order the data frame by the ranking order
  
  return(data)
} # Step 0

# Takes the dataframe and the gene set
# Returns the updated dataframe and the sum of the overlapping genes
assessing_overlaping <- function(data, candidate_list){
  
  indicator_function <- c() 
  sum_in_set = 0 # Keep track of the total rank metric for the genes overlapping
  
  for(i in data$ranking){ 
    if(data$names[i] %in% candidate_list){
      indicator_function <- c(indicator_function, 1) 
      sum_in_set = sum_in_set + abs(data$values[i])
    } else{
      indicator_function <- c(indicator_function, 0) 
    }
  }
  
  data["indicator"] <- indicator_function # Add vector as a new column of the dataframe
  
  # Get the counts of overlapping genes
  # hit <- sum(with(data, indicator == 1)) 
  # misses <- sum(with(data, indicator == 0)) 
  
  # Create a list to allow multiple returns
  my_list <- list("data" = data, "sum_in_set" = sum_in_set)
  
  return(my_list)
} # Step 1

# Takes the dataframe and the sum of the overlapping genes
# Returns the updated dataframe and the MES rank 
computing_sk_ES <- function(data, sum_in_set){
  
  misses <- sum(with(data, indicator == 0))
  
  running_ES <- c()
  
  ES = 0
  MES_value = 0 
  
  for(i in data$ranking){ 
    if(data$indicator[i] == 1){
      ES = ES + (abs(data$values[i]) / sum_in_set)
      running_ES <- c(running_ES, ES)
    } else {
      ES = ES - (1 / misses)
      running_ES <- c(running_ES, ES)
    }
    if(abs(ES) > abs(MES_value)){ 
      MES_value = ES
    }
  }
  
  data["running_ES"] <- running_ES 
  
  MES_rank <- which(data$running_ES == MES_value) # To know the ranking of the MES
  
  my_list <- list("data" = data, "MES_rank" = MES_rank, "MES_value" = MES_value)
  
  return(my_list)
} # Step 2 and 3

# Takes as input the dataframe and the number of permutations to be completed
# Returns a vector with the MES of each permutation
phenotype_permutation <- function(data_permutation, permutation){
  
  # Create a copy of the dataframe that contains the names, values and the indicator
  data_permuted <- data_permutation[,1:2]
  data_permuted['indicator'] <- data_permutation[,4]
  
  ES_null <- c() # Vector to keep the MES computed in each permutation
  pos_null <- c()
  neg_null <- c()
  
  for(i in 1:permutation){ # (5.3)
    
    # We make the permutation over the values
    # In this way we do not need to assess the overlaping
    data_permuted$values <- sample(data_permuted$values) # The sample() function allows to perform random permutations 
    
    # We need to get the absolute sum of the rank metric for those genes that overlap
    permuted_sum_in_set <- sum(abs(data_permuted[which(data_permuted$indicator == 1), 2]))
    
    # We directly make the ranking of the new dataframe
    data_permuted <- data_permuted[order(data_permuted$values),] # Ordering the data frame by the value of expression
    data_permuted["ranking"] <- c(length(data_permuted$names):1)
    data_permuted <- data_permuted[order(data_permuted$ranking),] # Order the data frame by the ranking order
    
    misses <- sum(with(data_permuted, indicator == 0))
    
    ES = 0
    MES_value = 0 
    
    for(i in data_permuted$ranking){ 
      if(data_permuted$indicator[i] == 1){
        ES = ES + abs(data_permuted$values[i]) / permuted_sum_in_set
      } else {
        ES = ES - (1 / misses)
      }
      if(abs(ES) > abs(MES_value)){ 
        MES_value = ES
      }
    }
    
    ES_null <- c(ES_null, MES_value)
    
    if(MES_value > 0){
      pos_null <- c(pos_null, MES_value)
    } else{
      neg_null <- c(neg_null, MES_value)
    }
  }
  
  mylist <- list("ES_null" = ES_null, "pos_null" = pos_null, "neg_null" = neg_null)
  
  return(mylist)
} # Step 5.1

# Main function call
# Returns the final dataframe, the MES rank (order) and the MES value (score)
OurGSEA_plus <- function(data, current_candidate){
  
  overlap <- assessing_overlaping(data, current_candidate) # Step 1
  
  data_sk_ES <- computing_sk_ES(overlap$data, overlap$sum_in_set) # Step 2 and 3
  
  mylist <- list("data" = data_sk_ES$data, "MES_rank" = data_sk_ES$MES_rank, "MES_value" = data_sk_ES$MES_value)
  
  return(mylist)
}

##########################################################
##########################################################
##########################################################

# Loading necessary data

###
# HALLMARK COLLECTION FOR HUMAN
###

# Loading the file
hallmark_data_human <- read.csv("h.all.v7.5.1.symbols.gmt", sep = "", header = FALSE)

transposed_data_human <- t(hallmark_data_human) # Change orientation of the original table (put rows as columns, and the other way round)
transposed_data_human <- as.data.frame(transposed_data_human)# Put the previous output as a dataframe. This table will have the categories as columns and the genes as rows

hallmark_human <- c()
multiple_description <- c()

for(i in 1:ncol(transposed_data_human)){
  hallmark_human <- c(hallmark_human, list(transposed_data_human[3:nrow(transposed_data_human),i]))
  multiple_description <- c(multiple_description, transposed_data_human[1,i])
}

###
# HALLMARK COLLECTION FOR MOUSE (translated from the human one)
###

# Load the translator table
# 2 columns file that contains genes correspondence between human and mouse
translator <- read.table("mart_export.txt", sep = ",", header = TRUE)

# Remove those rows that do not have a translator in mouse
filtered_translator <- translator[-which(translator$Gene.name == ""),]
filtered_translator <- filtered_translator[-which(filtered_translator$Human.gene.name == ""),]

hallmark_mouse <- c() # Generate the hallmark for mouse based on the translator table

for(cat in hallmark_human){ # Iterate through each gene set of the hallmark
  current_list <- c() # Vector to keep the translated gene of each set (empty for each new iteration)
  for(gene in cat){ # Iterate through each gene inside each gene set
    if(gene %in% filtered_translator$Human.gene.name){ # Check if the gene in the 
      candidates <- c(which(filtered_translator$Human.gene.name == gene))
      current_list <- c(current_list, filtered_translator$Gene.name[candidates[1]])
    }
    # The genes that do not have a correspondence are not included
  }
  candidate_mouse <- list(unique(current_list)) # Erase the genes that are repeated and put them together to a list
  hallmark_mouse <- c(hallmark_mouse, candidate_mouse) # Create a list of list, so that each gene set contains each translated genes
}
