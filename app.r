#general part, run only at app start


library(shiny)
library(shinydashboard)
library(gplots)
library(RColorBrewer)
library(ggplot2)
library(reshape)
library(ComplexHeatmap)
library(circlize)
library(viridis)	#for color palettes
library(rje)		#for color palettes
library(gridExtra)
library(plotly)
library(getopt) #for complexheatmap_single.r



source("helpers.R")
#source("complexheatmap_single.r")


#General part, done once
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#Data Input---------------------------------------------------------------------------
#Load data


table1 <- read.delim("data/normed_counts_orderd_development_ZB_Sven3.tsv", header=TRUE, check.names=F)
#table1 <- read.delim("data/matrix_log2fc.txt", header = TRUE, check.names = FALSE)
#table1=read.delim(file="data/normed_counts_orderd_development_ZB_Sven3.tsv",sep="\t",header=T,stringsAsFactors=T,row.names=1,check.names=FALSE)		#with header, column 0 = rownames, do not convert strings to character vectors
#reps=read.delim(file="data/reps.txt",sep="\t",header=F,stringsAsFactors=T,row.names=NULL,check.names=FALSE)		#with header, column 0 = rownames, do not convert strings to character vectors

#dynamic reps
reps <- data.frame(V1 = sub("(.*)_\\d$", "\\1", colnames(table1)[3:ncol(table1)]), V2 = colnames(table1)[3:ncol(table1)])

row.names(table1) <- table1[,1]



#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#UI---------------------------------------------------------------------------

ui <- dashboardPage(
  #Header---------------------------------------------------------------------------
  #Header Name upper left
 
   dashboardHeader(title = "MPI Result Viewer",
    
  #---------------------------------------------------------------------------
  #context list: notifications
  dropdownMenu(type = "notifications",
               
               notificationItem(
                 text = "Arrange ICON",
                 icon("warning")
                ),
               
               notificationItem(
                 text = "12 items delivered",
                 icon("truck"),
                 status = "success"
                )
              ),
  
  dropdownMenu(type = "tasks", badgeStatus = "success",
               
               taskItem(value = 5, color = "green",
                        "Documentation"
                ),
               taskItem(value = 17, color = "aqua",
                        "Heatmap"
                ),
               taskItem(value = 10, color = "yellow",
                        "enrichment"
                ),
               taskItem(value = 5, color = "red",
                        "Overall project"
                )
              )
  
            ),
  
  
  
  
  dashboardSidebar(
    #Side Bar-----------------------------------------------------------------------
    #Side Bar
    sidebarMenu(
      div(img(src="icon.png", height = 70, width = 70),em("Bioinformatics Core Unit")),
      
        menuItem("Overview", tabName = "overview", icon = icon("dashboard")),
        menuItem("Scatters", tabName = "scatter", icon = icon("area-chart")),
        menuItem("Heatmap", tabName = "heatmap", icon = icon("th"), selected = TRUE), # selected needs to be removed
        menuItem("Geneview", tabName = "genview", icon = icon("bar-chart")),
        menuItem("Enrichment", tabName = "enrichment", icon = icon("cc-mastercard"))
      
    )
  ),
  
  
  body<-dashboardBody(tabItems(
    
    #Body-----------------------------------------------------------------------
    #Main page for Overview
    tabItem(tabName = "overview",
            fluidRow(
              box(width=11,
                h2("General informations"), 
                #img(src="icon.png", height = 40, width = 40)
                p("hier kommt der fliesstext")
                ),
              box(width=1,
                  p(img(src="icon.png", height = 50, width = 50),
                  em("Bioinformatics Core Unit"))
              ),
              
              box(
                h2("noch a box")
                )
              
            )
    ),
    
    #scatter-----------------------------------------------------------------------
    #Main page for scatter
    tabItem(tabName = "scatter",
            fluidRow(
              box(width=12,title="Inputs",id="scatter_inputs",
                  column(4,
                         selectInput("scatter_xaxis",label="X axis:", choices = names(table1), selected="zygote_1")
                  ),
                  column(4,
                         selectInput("scatter_yaxis",label="Y axis:",choices = names(table1), selected="zygote_2")
                  )
              )
            ),
            
            fluidRow(
              box(width=12, title="Output", id="Output_scatter",
                plotOutput("plot_scatter",height="100%")
              )
            )
    ),
    #heatmap-----------------------------------------------------------------------
    #Main page for heatmap
    tabItem(tabName = "heatmap",
            fluidRow(
              box(width=12,title="Inputs",id="heatmap_inputs",
                  column(2,
                         selectInput("heat_select",label="Genes:",multiple=TRUE, choices = as.list(sort(unique(table1[,2]))), selected=table1[1,2]),
                    br(),
                         selectInput("heat_mode",label="Data transformation:",c("raw","log2","zscore"),selected="condition")
                  ),
                  column(2,
                         selectInput("heat_distrib",label="Data distribution",c("auto","one-sided","two-sided"), selected="auto"),
                    br(),
                         selectInput("heat_color",label="Color Type:",c("reds","viridis","plasma","inferno","magma", "blues", "heat", "cubehelix", "ylgnbu"),selected="BrBG")
                  ),
                  column(2,
                         selectInput("heat_clustering",label="Clustering",c("none", "row", "column", "both"), selected="both"), #both
                         br(),
                         selectInput("heat_clusterdist",label="Cluster Distance",c("euclidean", "pearson", "spearman", "kendall", "maximum", "manhattan", "canberra", "binary", "minkowski"), selected="manhattan")
                  ),
                  column(2,
                         selectInput("heat_clustermethod",label="Cluster method",c("average", "ward.D", "ward.D2", "single", "complete", "mcquitty", "median", "centroid"), selected="average"),
                         br(),
                         textInput("heat_unitlabel",label="Unit label",value = "Enter unit...")
                  ),
                  column(2,
                         checkboxInput("heat_reverse",label="Reverse Coloring:",value = FALSE),
                         br(),
                         checkboxInput("heat_rowlabel",label="Row Label",value = T) #true
                  ),
                  column(2,
                         checkboxInput("heat_columnlabel",label="Column Label",value = T) #true
                  )
              )
            ),
            
          
            #-u|--unitlabel       Label for unit displayed above color key.
            
            
            fluidRow(
              tabBox(width=12, title="Output", id="Output_heatmap",side="right",
                     tabPanel("Plot", plotOutput("plot_heatmap",height="100%")),
                     tabPanel("Table",dataTableOutput("table_heatmap"))
              )
            )
    ),
    #genview-----------------------------------------------------------------------
    ##Main page for genview
    #Layout: 2 rows, upper row with inputs, lower row with resulting graphs
    tabItem(tabName = "genview",
            fluidRow(
              box(width=12,title="Inputs",id="genview_inputs",
                  column(2,
                  selectInput("gv_select",label="Genes:",multiple=TRUE, choices = as.list(sort(unique(table1[,2]))), selected= table1[1,2])
                  ),
                  column(2,
                  selectInput("gv_facet",label="Grouping:",c("condition","gene"),selected="condition")
                  ),
                  column(2,
                  selectInput("gv_plottype",label="Type of Plot:",c("box","line","violin","bar"),selected="line")
                  ),
                  column(2,
                  selectInput("gv_color",label="Color Type:",c("Dark2","Accent","Reds","BrBG"),selected="BrBG")
                  ),
                  column(2,
                  sliderInput("gv_cols",label="Plot Columns",min = 1, max = 7, value =2)
                  ),
                  column(2,
                  sliderInput("gv_reserve",label="Reserve Input",min = 1, max = 7, value =2)
                  )
              )
            ),
            fluidRow(
            tabBox(width=12, title="Output", id="Output_genview",side="right",
                 #tabPanel("Plot", plotOutput("plot_genview",height="100%")),
                 tabPanel("Plot", plotlyOutput("plot_genview",height="100%")),
                 tabPanel("Table",dataTableOutput("table_genview"))
              )
            )
            
            
    ),
    #enrichment-----------------------------------------------------------------------
    #Main page for enrichment
    tabItem(tabName = "enrichment",
            h2("Gene sets and enrichment"),
            img(src="icon.png")
    )
    #-----------------------------------------------------------------------
    
  )
  )
)



# Server --------------------------------------------------------------------

server <- function(input, output, session) {
   #genview-------------------------------------------------------------
  #Section genview
  
  output$table_genview<-renderDataTable({
    genes_t<-table1[table1[[2]] %in% input$gv_select,1]
    values_t<-as.data.frame(table1[genes_t,])
    values_t
  },options=list(scrollX=TRUE))
  
  #output$plot_genview<-renderPlot({
  #  genes<-table1[table1$symbol %in% input$gv_select,1]
  #  values<-as.data.frame(table1[genes,])
  #  values2<-values[,3:92]
  #  width=1;
    #Function Call for Genview Plots
  #  dynamic_matrixsplit(values2,reps, input$gv_plottype, input$gv_facet,input$gv_color,input$gv_cols, width,input$gv_height)
  #  
  #}, height=900)

    
  output$plot_genview<-renderPlotly({
    genes<-table1[table1[[2]] %in% input$gv_select,1]
    values<-as.data.frame(table1[genes,])
    values2<-values[,3:ncol(table1)]
    width=1;
    #Function Call for Genview Plots
    x<-dynamic_matrixsplit(values2,reps, input$gv_plottype, input$gv_facet,input$gv_color,input$gv_cols, width,input$gv_height)
    ggplotly(x,height=900)
  })#, height=900)
  
  #heatmap-------------------------------------------------------------
  #Section heatmap
  source("helpers.R")  #removed in future
  output$table_heatmap<-renderDataTable({
    genes_t<-table1[table1[[2]] %in% input$heat_select,1]
    values_t<-as.data.frame(table1[genes_t,])
  },options=list(scrollX=TRUE))
  
  output$plot_heatmap<-renderPlot({
    genes<-table1[table1[[2]] %in% input$heat_select,1]
    values<-as.data.frame(table1[genes,])
    values_Heat<-values[,3:ncol(table1)]
    width=1
    
    
    #Function Call for Heatmap Plots
    #color_vector != null
    create_complexheatmap(m=values_Heat,mode=input$heat_mode, unitlabel=input$heat_unitlabel, rowlabel=input$heat_rowlabel, collabel=input$heat_columnlabel, clustering=input$heat_clustering, clustdist=input$heat_clusterdist, clustmethod=input$heat_clustermethod, distribution=input$heat_distrib, color_vector_onesided=input$heat_color, color_vector_twosided=input$heat_color, reverse_coloring=input$heat_reverse, optimize=T)
    
    
    ##single file | not working
    #create_complexheatmap(m, mode = input$heat_mode, unitlabel = input$heat_unitlabel, rowlabel = input$heat_rowlabel, collabel = input$heat_columnlabel, clustering = input$heat_clustering, clustdist = input$heat_clusterdist)
  }, height = 500 )
  
  #switch colors one-/two-sided
  observe({
    if (input$heat_distrib == "auto" & input$heat_mode != "raw" | input$heat_distrib == "two-sided") {
      updateSelectInput(session = session, inputId = "heat_color", choices = c("buwtrd", "rdblgr", "ylwtpu", "spectral"), selected = "spectral")
      
      #change label
      if(input$heat_mode == "zscore"){
        updateTextInput(session = session, inputId = "heat_unitlabel", value = "zscore")
      }else if(input$heat_mode == "log2"){
        updateTextInput(session = session, inputId = "heat_unitlabel", value = "log2")
      }
      
    }else{
      updateSelectInput(session = session, inputId = "heat_color", choices = c("reds","viridis","plasma","inferno","magma", "blues", "heat", "cubehelix", "ylgnbu"), selected = "reds")
      
      #change label
      if(input$heat_mode == "raw"){
        updateTextInput(session = session, inputId = "heat_unitlabel", value = "Enter unit...")
      }
    }
  
  })
  
    
  #scatter-------------------------------------------------------------
  #Section scatter
  
  output$plot_scatter<-renderPlot({
    selectedData <- table1[, c(input$scatter_xaxis, input$scatter_yaxis)]
    write.csv(selectedData,file="test.txt")
      
    
    ggplot(table1,aes(x=input$scatter_xaxis,y=input$scatter_yaxis)) + geom_point() 
    
  }, height=400)
}

shinyApp(ui, server)