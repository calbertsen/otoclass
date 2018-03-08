
library(shiny)
library(rmarkdown)
render("method.Rmd","html_document",runtime="shiny")


##UI
shinyUI(
    fluidPage(titlePanel("Otolith Classification: Bias correction in proportion estimates"),
              sidebarLayout(
                  sidebarPanel(
                      sliderInput("Nb",
                                  "Number of bootstrap simulations:",
                                  min = 1,
                                  max = 500,
                                  value = 5,
                                  step = 1,
                            round = FALSE),
                      sliderInput("Ns",
                                  "Number of simulations:",
                                  min = 100,
                                  max = 1000,
                                  value = 400,
                                  step = 1,
                                  round = FALSE),
                      sliderInput("ps",
                                  "Proportion in first group for training data:",
                                  min = 0.05,
                                  max = 0.95,
                                  value = 0.5,
                                  step = 0.01,
                                  round = FALSE),
                      sliderInput("pstest",
                                  "Proportion in first group for test data:",
                                  min = 0.05,
                                  max = 0.95,
                                  value = 0.5,
                                  step = 0.01,
                                  round = FALSE),
                      sliderInput("x1",
                                  "X mean for first group:",
                                  min = -10,
                                  max = 10,
                                  value = -1,
                            step = 0.1,
                                  round = FALSE),
                      sliderInput("y1",
                                  "Y mean for first group:",
                                  min = -10,
                                  max = 10,
                                  value = -1,
                                  step = 0.1,
                                  round = FALSE),
                      sliderInput("sd",
                                  "Scale coefficient for first group relative to second group:",
                                  min = 0.1,
                                  max = 5,
                                  value = 1,
                                  step = 0.1,
                                  round = FALSE),
                      sliderInput("rho",
                                  "Correlation between covariates:",
                                  min = -1,
                                  max = 1,
                                  value = 0,
                                  step = 0.01,
                                  round = FALSE)
                  ),
                  mainPanel(
                      tabsetPanel(
                          tabPanel("Data",
                                   plotOutput("dataplot")
                                   ),
                          tabPanel("Bootstrap Plots",
                                   h2("Estimated proportions with correction"),
                                   plotOutput("proportionplot"),
                                   h2("Estimated proportions before correction"),
                                   plotOutput("rawproportionplot"),
                                   h2("Estimated confusion matrix"),
                                   plotOutput("confusionmat")
                                   ),
                          tabPanel("Bootstrap Summary",
                                   h2("Summary of estimates with correction"),
                                   tableOutput("summarytable"),
                                   h2("Summary of estimates before correction"),
                                   tableOutput("rawsummarytable")
                                   ),
                          tabPanel("Method",
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

