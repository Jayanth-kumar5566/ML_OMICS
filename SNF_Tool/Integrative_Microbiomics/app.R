library(shiny)
library(SNFtool)
library("vegan")
source("functions.R")
ui <- fluidPage(
  headerPanel("Integrative Microbiomics",windowTitle="Microbiomics"),
  navbarPage("Integrative Microbiomics",
             tabPanel("Biome Submission",fluidRow(
               column(4,fileInput(inputId = "biome1",label = "Biome: 1"),
                      plotOutput(outputId = "biome1_plot"),
                      numericInput("k_biome1", label = h3("Number of Clusters"), value = 2),
                      actionButton(inputId = "go_biome1",label = "Upload Biome: 1")),
               column(4,fileInput(inputId = "biome2",label = "Biome: 2"),
                      plotOutput(outputId = "biome2_plot"),
                      numericInput("k_biome2", label = h3("Number of Clusters"), value = 2),
                      actionButton(inputId = "go_biome2",label = "Upload Biome: 2")),
               column(4,fileInput(inputId = "biome3",label = "Biome: 3"),
                      plotOutput(outputId = "biome3_plot"),
                      numericInput("k_biome3", label = h3("Number of Clusters"), value = 2),
                      actionButton(inputId = "go_biome3",label = "Upload Biome: 3")))
             ),
             tabPanel("Parameter Selection",
                      actionButton(inputId = "go",label = "Merge")),
             tabPanel("Cluster Visuvalization",
                      plotOutput(outputId = "Cluster_matrix"),
                      downloadButton(outputId = "download",label="Downlod Label Files")))
)

server <- function(input, output) {
  data1=eventReactive(input$go_biome1,{read.csv(input$biome1$datapath,header = TRUE,row.names = 1)})
  output$biome1_plot<-renderPlot({biome_plot(data1(),input$k_biome1)})
  data2=eventReactive(input$go_biome2,{read.csv(input$biome2$datapath,header = TRUE,row.names = 1)})
  output$biome2_plot<-renderPlot({biome_plot(data2(),input$k_biome2)})
  data3=eventReactive(input$go_biome3,{read.csv(input$biome3$datapath,header = TRUE,row.names = 1)})
  output$biome3_plot<-renderPlot({biome_plot(data3(),input$k_biome3)})
  
  # data2=eventReactive(input$go_biome2,{read.csv(input$biome2$datapath)})
  # data3=eventReactive(input$go_biome3,{read.csv(input$biome3$datapath)})
  # W=eventReactive(input$go,{
  #   d1_dsim=vegdist(data1()[,-1],method='bray',diag=TRUE,upper=TRUE) #-1 to remove index
  #   d2_dsim=vegdist(data2()[,-1],method='bray',diag=TRUE,upper=TRUE) #-1 to remove index
  #   W1=(as.matrix(d1_dsim)-1)*-1
  #   W2=(as.matrix(d2_dsim)-1)*-1
  #   W = SNF(list(W1,W2),23,15)
  # })
  # label=reactive(spectralClustering(W(), K = 2))
  # output$Cluster_matrix<-renderPlot({
  #   displayClusters(W(),label())
  # })
  # output$download<-downloadHandler(
  #   filename = function() {
  #     paste('labels-', Sys.Date(), '.csv', sep='')
  #   },
  #   content = function(con) {
  #     write.csv(label(), con)
  #   }
  # )
}
shinyApp(ui = ui, server = server)