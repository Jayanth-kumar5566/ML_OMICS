library(shiny)
library(SNFtool)
library("vegan")
ui <- fluidPage(
  headerPanel("Integrative Microbiomics",windowTitle="Microbiomics"),
  fileInput(inputId = "biome1",label = "Biome: 1"),
  fileInput(inputId = "biome2",label = "Biome: 2"),
  actionButton(inputId = "go",label = "Submit the Biomes"),
  plotOutput(outputId = "Cluster_matrix"),
  downloadButton(outputId = "download",label="Downlod Label Files")
)

server <- function(input, output) {
  data1=eventReactive(input$go,{read.csv(input$biome1$datapath)})
  data2=eventReactive(input$go,{read.csv(input$biome2$datapath)})
  W=eventReactive(input$go,{
    d1_dsim=vegdist(data1()[,-1],method='bray',diag=TRUE,upper=TRUE) #-1 to remove index
    d2_dsim=vegdist(data2()[,-1],method='bray',diag=TRUE,upper=TRUE) #-1 to remove index
    W1=(as.matrix(d1_dsim)-1)*-1
    W2=(as.matrix(d2_dsim)-1)*-1
    W = SNF(list(W1,W2),23,15)
  })
  label=reactive(spectralClustering(W(), K = 2))
  output$Cluster_matrix<-renderPlot({
    displayClusters(W(),label())
  })
  output$download<-downloadHandler(
    filename = function() {
      paste('labels-', Sys.Date(), '.csv', sep='')
    },
    content = function(con) {
      write.csv(label(), con)
    }
  )
}
shinyApp(ui = ui, server = server)