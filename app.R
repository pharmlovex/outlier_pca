#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
source("helper.R")
library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Data Cleaning With PCA"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            fileInput("datafile",
                      "Upload count data:",
                      accept = c(
                        "text/csv",
                        "test/comma-separated-values,text/plain",
                        ".tsv",
                        ".csv")
            )
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("scaleplot"),
           tableOutput("gene_wt"),
           tableOutput("sample_pc"),
           plotOutput("variance_plot"),
           plotOutput("pca_plot"),
           plotOutput("outlier_plot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  result <- reactive({
      infile <- input$datafile
      if (is.null(infile))
        return(NULL)
      wrangle(infile$datapath)
  })
  output$scaleplot <- renderPlot({
        result = result()
        grid.arrange(result$dp_before_scale, result$dp_after_scale, ncol =2)
    })
  output$gene_wt <- renderTable({
    result = result()
    result$gene_pc_df
  })
  output$sample_pc <- renderTable({
    result = result()
    result$sample_pc_df
  })
  output$variance_plot <- renderPlotly({
    result = result()
    ggplotly(result$variance_explained_plot)
  })
  output$pca_plot <- renderPlot({
    result = result()
    result$pca_plot
  })
  output$outlier_plot <- renderPlot({
    result = result()
    result$joint_prob_density_plot
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
