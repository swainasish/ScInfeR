library(shiny)
library(shinythemes)
library(ggplot2)
library(htmlwidgets)
library(htmltools)
library(networkD3)
library(dplyr)
library(readxl)
library(markdown)
library(plotly)
library(googlesheets4)
gs4_deauth()


all_tissu_types <- c("Bladder", "Bone marrow", "Brain", "Breast", "Eye", "Heart", 
                     "Intestine", "Kidney", "Liver", "Lungs", "Pancreas", "PBMC", "Skin")
figshare_url <- read.csv("datasets/download_url_csv.csv")
rownames(figshare_url) <- figshare_url$Tissue
url <- "https://docs.google.com/spreadsheets/d/1qOawAgIIESHAXkAkN8xCHH65m3T60ytDYaJgyzM2Wrs/edit?usp=sharing"
suggest_marker_df <- read_sheet(url)
colnames(suggest_marker_df) <- c("Time-stamp","User name","Tissue","Cell type","Marker name","Weights")

# marker_hs<- read.csv("scinfer_combined_hs.csv",
#                      row.names = 1)
# links<-read_excel("hirarchy_tissue_scinfer.xlsx",
#                      sheet = "Lungs")
# 
# subsetdf=read_excel("datasets/scinferdb_MASTER_UMAP.xlsx",
#                     sheet = "Bladder",n_max = 5000)
# row.names(subsetdf) <- NULL
# ggplot(data=subsetdf,aes(x=umap_1,y=umap_2,color=Celltype))+geom_point(size=0.2)
# From these flows we need to create a node data frame: it lists every entities involved in the flow
# nodes <- data.frame(
#   name=c(as.character(links$source), 
#          as.character(links$target)) %>% unique()
# )
# 
# # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
# links$IDsource <- match(links$source, nodes$name)-1 
# links$IDtarget <- match(links$target, nodes$name)-1

# Make the Network
# p <- sankeyNetwork(Links = links, Nodes = nodes,
#                    Source = "IDsource", Target = "IDtarget",
#                    Value = "value", NodeID = "name", 
#                    sinksRight=FALSE)

CSS <- "
.button {
  font: bold 11px Arial;
  text-decoration: none;
  background-color: #EEEEEE;
  color: #333333;
  padding: 4px 8px 4px 8px;
  border-top: 1px solid #CCCCCC;
  border-right: 1px solid #333333;
  border-bottom: 1px solid #333333;
  border-left: 1px solid #CCCCCC;
}"



ui <- fluidPage(theme = shinytheme("yeti"),
                tags$head(
                  tags$style(HTML(CSS)) ),
                navbarPage(
                  # theme = "cerulean",  # <--- To use a theme, uncomment this
                  "ScInfeR",
                   navbarMenu("Documentation",
                              tabPanel("Installation",
                                       fluidRow(
                                         htmltools::includeMarkdown("about.md") )),
                              tabPanel("Annotate scRNA-seq (Marker-based)",
                                       fluidRow(
                                         htmltools::includeMarkdown("about.md") )),
                              tabPanel("Annotate scRNA-seq (Reference-based)",
                                       fluidRow(
                                         htmltools::includeMarkdown("about.md") )),
                              tabPanel("Annotate scRNA-seq (Hybrid-based) [Recommended]",
                                       fluidRow(
                                         htmltools::includeMarkdown("about.md") )),
                              tabPanel("Installation",
                                       fluidRow(
                                         htmltools::includeMarkdown("about.md") )),
                              tabPanel("Installation",
                                       fluidRow(
                                         htmltools::includeMarkdown("about.md") )),
                              tabPanel("Installation",
                                       fluidRow(
                                         htmltools::includeMarkdown("about.md") )),
                              tabPanel("Installation",
                                       fluidRow(
                                         htmltools::includeMarkdown("about.md") )),
                              tabPanel("Installation",
                                       fluidRow(
                                         htmltools::includeMarkdown("about.md") )),
                              tabPanel("Installation",
                                       fluidRow(
                                         htmltools::includeMarkdown("about.md") ))),
                  
                  tabPanel("CellMarkerDB",
                           sidebarPanel(
                             selectizeInput("tissue_type", "Choose tissue type", 
                                            choices = all_tissu_types),
                             sankeyNetworkOutput("hirar_plot")),
                           mainPanel(
                             h3(textOutput("txt_out1")),
                             fluidRow(
                               splitLayout(
                                           DT::dataTableOutput("table")),
                               downloadButton("downloadData", "Download")
                               ))),
                  tabPanel("scRNA-seqDB",
                           sidebarPanel(
                             selectizeInput("tissue_type_ref", "Choose tissue type",
                                            choices = all_tissu_types))
                           ,mainPanel(h6("This UMAP is a representation dataset from the full size scRNA-seq reference"),
                             plotlyOutput("umap_plot",height = "600px", width = "700px"),
      
                                      downloadButton("downloadData_scref", "Download (scRNA-seq)"))),
                  tabPanel("Suggest a marker",
                           sidebarPanel(htmlOutput("map")),
                           mainPanel(h3("Already suggested markers by the users"),
                             DT::dataTableOutput("suggest_dt"))),
                  tabPanel("Citation",
                           "Blank")
                            ))

# Define server function  
server <- function(input, output) {
  datasetInput <- reactive({
    switch(input$tissue_type,
           cardf[cardf$cyl==input$tissue_type,])
  })
  
  subsetdf <- reactive({subsetdf=read_excel("datasets/scinfer_combined_hs.xlsx",
                                            sheet = input$tissue_type)
  row.names(subsetdf) <- NULL
  subsetdf
  })
  # Table of selected dataset ----
  output$table <- DT::renderDataTable(DT::datatable({
    subsetdf=read_excel("datasets/scinfer_combined_hs.xlsx",
                        sheet = input$tissue_type)
    row.names(subsetdf) <- NULL
    subsetdf
  }))
  output$suggest_dt <- DT::renderDataTable(DT::datatable({
    suggest_marker_df
  }))
  

output$umap_plot <- renderPlotly({
  umapdf=read_excel("datasets/scinferdb_MASTER_UMAP.xlsx",
                      sheet = input$tissue_type_ref,n_max=8000)
  row.names(umapdf) <- NULL
  umapplt <- ggplot(data=umapdf,aes(x=umap_1,y=umap_2,color=Celltype))+geom_point(size=0.4)+xlab("UMAP1")+ylab("UMAP2")+theme_bw()
  umapplt})
  
output$hirar_plot <- renderSankeyNetwork({
    links<-read_excel("datasets/hirarchy_tissue_scinfer.xlsx",
                      sheet = input$tissue_type)
    nodes<- data.frame(
      name=c(as.character(links$source), 
             as.character(links$target)) %>% unique())
    
    links$IDsource <- match(links$source, nodes$name)-1 
    links$IDtarget <- match(links$target, nodes$name)-1 
    splot<- sankeyNetwork(Links = links, Nodes = nodes,
                  Source = "IDsource", Target = "IDtarget",
                  Value = "value", NodeID = "name", 
                  sinksRight=FALSE,fontSize = 12)
    splot
  })
  
  subsetdf <- reactive({subsetdf=read_excel("datasets/scinfer_combined_hs.xlsx",
                                            sheet = input$tissue_type)
  row.names(subsetdf) <- NULL
  subsetdf
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("ScInfeR","_",input$tissue_type,'.csv')
    },
    content = function(file) {
      write.csv(subsetdf(), file, row.names = FALSE)
    }
  )
  output$downloadData_scref <- downloadHandler(
    filename = function() {
      paste0("ScInfeR","_",input$tissue_type_ref,'.rds')
    },
    content = function(file) {
      url_name = figshare_url[input$tissue_type_ref,"url"]
      browseURL(url_name)
    }
  )
  #output$downloadData_scrna <- downloadLink("https://figshare.com/ndownloader/files/34702051")
  
  output$txt_out1<-renderText({paste("Cell marker list for tissue:",input$tissue_type)})
  output$txtout <- renderText({
    paste( input$txt1, input$txt2, sep = " " )
  })
  output$map <- renderUI({
    tags$iframe(seamless="seamless",
                src= "https://docs.google.com/forms/d/e/1FAIpQLSf0ApBMYzQWez5DFZI_a-DoGDYvKvUGyHoh6Ox9JOEJ4B726Q/viewform?embedded=true",
                width=350,
                height=800)
  })
} # server



# Create Shiny object
shinyApp(ui = ui, server = server)
