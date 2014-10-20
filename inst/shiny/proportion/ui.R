
library(shiny)
library(rmarkdown)
render("method.Rmd","html_document",runtime="shiny")


##UI
shinyUI(
    fluidPage(titlePanel("Otolith Classification: Estimating proportions"),
              sidebarLayout(
                  sidebarPanel(
                      helpText("If you do not have a csv file to upload, you can download an example and upload that:"),
                      downloadButton('downloadData', 'Download'),
                      fileInput('datfile', 'Choose file to upload (Please upload a file with as few columns as possible; max 60 MB):',
                                accept = c(
                                    'text/csv',
                                    'text/comma-separated-values',
                                    'text/tab-separated-values',
                                    'text/plain',
                                    '.csv',
                                    '.tsv'
                                )
                                ),
                      checkboxInput('header', 'Header', TRUE),
                      radioButtons('sep', 'Separator',
                                   c(Comma=',',
                                     Semicolon=';',
                                     Tab='\t'),
                                   ','),
                      radioButtons('quote', 'Quote:',
                                   c(None='',
                                     'Double Quote'='"',
                                     'Single Quote'="'"),
                                   '"'),
                      uiOutput("selectGroup"),
                      uiOutput("selectVariables"),
                      h4("Bootstrap options"),
                      sliderInput("Nb",
                                  "Number of bootstrap simulations:",
                                  min = 1,
                                  max = 500,
                                  value = 5,
                                  step = 1,
                                  round = FALSE),
                      h5("Number of individuals in training sample:"),
                      uiOutput("slidersTrain"),
                      h5("Number of individuals in test sample:"),
                      uiOutput("slidersTest") 
                  ),
                  mainPanel(
                      tabsetPanel(
                          tabPanel("Data",
                                   h2("Uploaded data"),
                                   tableOutput("uploadHead"),
                                   h2("Number of individuals in each group"),
                                   tableOutput("groupTable"),
                                   h2("Selected covariates"),
                                   plotOutput("plotVars"),
                                   h2("Bootstrap proportions"),
                                   tableOutput("bootprop")
                                   ),
                          tabPanel("Bootstrap Results",
                                   h2("Bootstrap results"),
                                   uiOutput('selectShowGroup'),
                                   h3("Estimated proportions"),
                                   plotOutput("proportionplot"),
                                   h2("Summary of estimates with correction"),
                                   tableOutput("summarytable"),
                                   h2("Classification success rate"),
                                   plotOutput("ccplot2tal")
                                   ),
                          tabPanel("Description",
                                   includeHTML("method.html")
                                   )
                      )   
                  )
              ),
              tags$div(style = "position: fixed; bottom: 0; text-align:center;width: 100%;",
                       p(list("Christoffer Moesgaard Albertsen, 2014, ",
                              a(href="mailto:cmoe@aqua.dtu.dk","cmoe@aqua.dtu.dk")
                              ))
                       )
              ))

