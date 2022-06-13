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

# load module functions
source("GSEAfunctions.R")
# load ui elements
source("ui.R")
# load server function
source("server.R")

# Run the application 
shinyApp(ui = ui, server = server)
# runApp(shinyApp(ui, server), launch.browser = TRUE)