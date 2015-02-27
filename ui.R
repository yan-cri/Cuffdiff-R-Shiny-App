library(shiny)

shinyUI(navbarPage("RNA-seq cuffdiff result visulization",
                 
      tabPanel("Data input",
               fluidPage(
                 fluidRow(
                   column(12, 
                          textInput(inputId = "file",label = h1("Data folder name"),value = "testdata"),
                           
                          helpText("Provide the folder name where you cuffdiff analysis results are saved inside the data directory, or the full path to the folder where your cuffdiff analysis results are saved")
                          ,
                          
                          actionButton("action", label = "Submit")
                          )
                   ),
                 fluidRow(
                           column(12,
                                  h2("Your cuffdiff analysis results are located at"),verbatimTextOutput("fileinput") 
                           ),
                           column(12,
                                  h3("cuffdiff analysis summary is as below:"),tableOutput("dataSummary")
                                  ),
                           column(12,
                                  h3("Sample summary is as below"),tableOutput("sampleSummary")
                                  )
                           
                           ),
                 br(),
                 br(),
                 br(),
                 mainPanel(position="bottom",img(src="cri.png"))
                 )
               ),
                   
      tabPanel("Sample Summary Plots",
               sidebarLayout(
                 sidebarPanel(
                   helpText("Click below button to chose whether use biological replicates"),
                   
                   selectInput("var", 
                               label = "Biological replicates or not",
                               choices = c("yes"=1, "no"=0),
                               selected = 1)
                 ),
                 
                 mainPanel(
                   tabsetPanel(
                     tabPanel("MDS plot", plotOutput("mdsPlot")),
                     tabPanel("Sample heatmap", plotOutput("sammpleHeatmap")),
                     tabPanel("Sample dispersion plot", plotOutput("dispersionPlot")),
                     tabPanel("Sample distribution plot", plotOutput("densPlot")),
                     tabPanel("Sample box-plot", plotOutput("sampleboxPlot")),
                     tabPanel("Sample dendergram", plotOutput("dendergramPlot")),
                     tabPanel("Sample pairwise scatterPlot", plotOutput("csScatterPlot"))
                   )
                 )
               )
        ),
      
      tabPanel("DEGs Summary",
               fluidPage(
                 fluidRow(column(4,
                                 h4("Fold Change (FC) and False Discovery Rate (FDR) cutoff option"),
                                 textInput(inputId = "degFC",label = "FC cutoff value",value = "1.5"), 
                                 textInput(inputId = "degFDR",label = "FDR cutoff value",value = "0.05"),
                                 br(),
                                 actionButton("degaction", label = "Submit")
                                 ),
                 column(4,
                        h4("No. DEGs with FC cutoff is as below:"), tableOutput("degSummary")
                 ),
                 column(4,
                        h4("List of DEGs names based on the provided cutoff criteria is as below:"), verbatimTextOutput("cuffdiffCommand"),
                        br()
                 )
                 ),
                 fluidPage(column(12,
                                  h3("DEGs Lists", align="center"),
                                  dataTableOutput("degList")
                                  )
                   )


                 )
               ),
      
      tabPanel("Gene of Interest",
               sidebarLayout(
                 sidebarPanel(
                   helpText("Please type in your genes of interest (GOI) list in the below box seperated by comma"),
                   
                   textInput(inputId = "goiList",label = "List of genes of interest",value = "FGFR1, RAB3A"),
                   
                   br(),
                   br(),
                   
                   helpText("Click below button to chose whether use biological replicates"),
                   
                   selectInput("goivar", 
                               label = "Biological replicates or not",
                               choices = c("yes"=1, "no"=0),
                               selected = 1)
                   ,
                   
                   actionButton("goiaction", label = "Submit")
                   
                   ),
                 
                 mainPanel(
                   tabsetPanel(
                     tabPanel("GOI heatmap", imageOutput("goiHeatmap")),
                     tabPanel("GOI overall expression", dataTableOutput("goiExpression")),
                     tabPanel("GOI summarized expression", br(), tableOutput("goiExpressionSummary")),
                     tabPanel("GOI boxplot", plotOutput("goiBoxPlot")),
                     tabPanel("GOI test results",  dataTableOutput("goiDiff")),
                     tabPanel("GOI FC barplot",  plotOutput("goiFCbarplot"))
                     )
                   )
                 )
               
        )
    

))

