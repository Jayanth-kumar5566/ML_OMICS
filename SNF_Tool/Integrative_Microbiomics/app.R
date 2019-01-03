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
             tabPanel("Parameter Selection",fluidRow(
               column(6,
                      h4("Choose a Similarity Measure/Metric"),
                      p("A metric/Similarity measure is a measure or function 
                        that defines distances or similarities between each sample.
                        For example, Bray-Curtis Similarity in 
                        case of ecological or Microbiomic samples ")),
               column(6,
                      h4("Choose a Merging Algorithm/Method"),
                      p("A Merging algorithm is a method through which the datasets
                        are integrated")
                      )),
               fluidRow(
                 hr(),
                 column(6,
                        selectInput("metric","Metric/Simialrity Measure :",c("Bray-Curtis")),
                        textOutput("metric_out"),
                        hr(),
                        actionButton(inputId = "merge",label = "Merge")),
                 column(6,
                        selectInput("method","Merging Method :",c("SNF")),
                        textOutput("method_out"),
                        uiOutput('ui'))
               )),
             tabPanel("Cluster Visuvalization",
                      fluidRow(
                        column(6,
                               h4("Modified Average Silhouette Width"),
                               p("The below table represents the average silhouette width for
                                 different number of clusters.The silhouette value is a measure
                                 of how similar an object is to its own cluster (cohesion)
                                 compared to other clusters (separation). The silhouette 
                                 ranges from −1 to +1, where a high value indicates that the
                                 object is well matched to its own cluster and poorly matched 
                                 to neighboring clusters. Since this is the average sihouette 
                                 width, relying soley on these values for number of cluster
                                 selection is not advisable.")),
                        column(6,
                               h4("Optimal Number of Clusters"),
                               p("Two different methods 1.Best Eigen Gap  2. Rotation cost
                                  are used  to identify the optimal/best number of clusters.
                                ")
                               )
                      ),
                      fluidRow(
                        column(6,tableOutput("k_Table")),
                        column(6,tableOutput("cluster_est"))
                      ),
                      fluidRow(
                        hr(),
                        column(6,
                               plotOutput(outputId = "merged_biome_plot")
                               ),
                        column(6,
                               numericInput("merged_biome", label = h4("Optimal Number of Clusters"),value=2),
                               downloadButton(outputId = "download",label="Downlod Label Files"))
                      )
                      )))
server <- function(input, output) {
  data1=eventReactive(input$go_biome1,{
    if (is.null(input$biome1$datapath))
    {return(NULL)}
    else{
      read.csv(input$biome1$datapath,header = TRUE,row.names = 1)
    }
  })
  output$biome1_plot<-renderPlot({biome_plot(data1(),input$k_biome1)})
  data2=eventReactive(input$go_biome2,{
    if (is.null(input$biome2$datapath))
    {return(NULL)}
    else{
      read.csv(input$biome2$datapath,header = TRUE,row.names = 1)
    }
  })
  output$biome2_plot<-renderPlot({biome_plot(data2(),input$k_biome2)})
  data3=eventReactive(input$go_biome3,{
    if (is.null(input$biome3$datapath))
    {return(NULL)}
    else{
      read.csv(input$biome3$datapath,header = TRUE,row.names = 1)
    }
  })
  output$biome3_plot<-renderPlot({biome_plot(data3(),input$k_biome3)})
  
output$metric_out<-renderText({
  if(input$metric=="Bray-Curtis"){"Bray-Cutis Similarity
  is a statistic used to quantify the compositional dissimilarity
  between two different sites"}})  
output$method_out<-renderText({
  if(input$method=="SNF"){"Similarity Network Fusion, 
    is a new computational method for data integration. 
    SNF first constructs a sample similarity network for each of 
    the data types and then iteratively integrates these networks using
    a novel network fusion method. Working in the sample network space 
    allows SNF to avoid dealing with different scale, collection bias and 
    noise in different data types. Integrating data in a non-linear fashion 
    allows SNF to take advantage of the common as well as complementary i
    nformation in different data types "}
  })

output$ui<-renderUI({
if(is.null(input$method))
  return()
  switch(input$method,
       "SNF"=tagList(
        numericInput("K_nn","K Neareast Neighbours",value = round(dim(data1())[1]/10)),
        numericInput("t_iter","Number of Iterations",value=20))
       )
  })

data_merge=eventReactive(input$merge,
                        {merge_snf(list(data1(),data2(),data3()),input$K_nn,input$t_iter)})

output$k_Table<-renderTable(max_k(data_merge()))
output$cluster_est<-renderTable({as.data.frame(estimateNumberOfClustersGivenGraph(data_merge()))})

output$merged_biome_plot<-renderPlot({biome_plot(data_merge(),input$merged_biome)})

label=reactive({label_create(data_merge(),input$merged_biome)})

output$download<-downloadHandler(
   filename = function() {
     paste('labels-', Sys.Date(), '.csv', sep='')
   },
   content = function(con) {
     write.csv(label(), con)
   },
   contentType = "csv"
 )



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
  
}
shinyApp(ui = ui, server = server)

#Todo
#Find out how to set the empty data3() to NULL
#Then change the merge_snf code appropriately
#Try rendering the table of silhouette and the plot 

