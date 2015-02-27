# Cuffdiff-R-Shiny-App

## Introduction
This git repo includes a R shiny App which can be used to visualize your obtained analysis results from cuffdiff saved in your local PC.

### Required software and packages to run this App
This is a R shiny App. In order to run this App, you need to have R and R "shiny" package installed. 

## Run this App in your local PC
1. Download this App to your local PC with gitclone or download zipped files.
2. Initial your R program.
3. Load shiny library in R by command: 
   library(shiny)
4. Run this App in R with the below comman:
   shiny::runApp("Path/to/download/Cuffdiff-R-Shiny-App/")
5. A html or a popup window (if use RStudio to run the App) will display the cuffdiff results from a test data this App initially is tested on.
6. Use the "Data folder name" in the Input navigation page to change the input data into your own cuffdiff results.
   Here, you need to provide the full path where you put your cuffdiff obtained results.
