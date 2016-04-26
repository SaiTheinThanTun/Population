#learning Reed-Forst
library(shiny)

ui <- fluidPage(
  fluidRow(h2("1-(1-x)^y"),
    column(4,
           h3("x value"),
           sliderInput(inputId = 'x',
                       label = 'x value',
                       value = .03, min = 0, max= 1)),
    column(4,
           h3('y values'),
           numericInput(inputId='y1',
                       label='y1 value',
                       value=2.5)),
    column(4,
           h3('y values'),
           numericInput(inputId='y2',
                        label='y2 value',
                        value=100))
  ),
  fluidRow(h2('output graph'),
           column(5,
                  plotOutput(outputId = 'graph')))
)

server <- function(input,output){
  output$graph <- renderPlot({
    f <- function(y){
      1-(1-input$x)^y
    }
    y <- seq(from=input$y1, to=input$y2, length=100)
    
    plot(y,f(y))
    
  })
}

shinyApp(ui = ui, server = server)
