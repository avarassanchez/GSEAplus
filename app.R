# Loading the required libraries

library(ggplot2)
library(shiny) 
library(shinydashboard)
library(knitr) 
library(data.table)
library(rintrojs) 
library(shinyWidgets)
library(rmarkdown) 

# library(dplyr)
# library(devtools)
# library(shinyFeedback)
# library(progress)
# library(cowplot) 

# load module functions
source("GSEAfunctions.R")
# load ui elements
source("ui.R")
# load server function
source("server.R")

# Run the application 
shinyApp(ui = ui, server = server)
# runApp(shinyApp(ui, server), launch.browser = TRUE)