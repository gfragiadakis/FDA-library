library(shiny)

shinyUI(fluidPage(

  # Application title
  titlePanel("Defining modules in adjacency matrices"),
  plotOutput(outputId = "plot", width = "100%", brush = brushOpts(id = "plot_brush", fill = "#ccc", direction = "x")),
  hr(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),br(),br(),br(),br(),br(),br(),br(),br(),

  fluidRow(

    column(
           verbatimTextOutput("info"), width = 6
    ),
    column(width = 3,
           h3("Documenting modules"),
           actionButton("save", label = "Save module")
    )
  ), br(),br(),br()
))


