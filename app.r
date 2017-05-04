#general part, run only at app start
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
library(data.table)
library(heatmaply)



source("helpers.R")
#source("complexheatmap_single.r")


#General part, done once
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#Data Input---------------------------------------------------------------------------
#Load data

table1 <- fread("data/normed_counts_orderd_development_ZB_Sven3_categorized.tsv", header = TRUE)
setkey(table1, id)

#print(head(table1))
# reps --------------------------------------------------------------------


#get categories
#reps <- reps_maker()

#columnnames unique?
if(any(duplicated(colnames(table1)))){
  #make columnnames unique
  allColumns <- colnames(table1)[sapply(table1, is.numeric)] #with duplicates
  table1 <- uniqueColumns(table1)
  
  reps <- data.table(allColumns, colnames(table1)[sapply(table1, is.numeric)])
  #print(reps)
}else{
  #dynamic reps
  reps <- data.table(V1 = sub("(.*)_\\d$", "\\1", colnames(table1)[sapply(table1, is.numeric)]), V2 = colnames(table1)[sapply(table1, is.numeric)])

}
setkey(reps, V2)


#table1 <- read.delim("data/normed_counts_orderd_development_ZB_Sven3_big.tsv", header=TRUE, check.names=F)
#table1 <- read.delim("data/matrix_log2fc.txt", header = TRUE, check.names = FALSE)
#table1=read.delim(file="data/normed_counts_orderd_development_ZB_Sven3.tsv",sep="\t",header=T,stringsAsFactors=T,row.names=1,check.names=FALSE)		#with header, column 0 = rownames, do not convert strings to character vectors
#reps=read.delim(file="data/reps.txt",sep="\t",header=F,stringsAsFactors=T,row.names=NULL,check.names=FALSE)		#with header, column 0 = rownames, do not convert strings to character vectors


# selectInputData ------------------------------------------------------------------------
genes <- sort(unique(table1[,2])[[1]])
columns_num <- sort(unique(colnames(table1)[sapply(table1, is.numeric)]))
columns <- sort(unique(colnames(table1)))





message("Data loaded")

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
        menuItem("Scatters", tabName = "scatter", icon = icon("area-chart"), 
                 menuSubItem(text = "Scatter", tabName = "scatter"),
                 menuSubItem(text = "Category", tabName = "scatter_cat")), # selected needs to be removed 
        menuItem("Heatmap", tabName = "heatmap", icon = icon("th"), selected = TRUE), 
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
                  column(3,
                         textInput(inputId = "scatter_X_label", label = "X axis label:", placeholder = "Custom label"),
                         selectInput("scatter_xaxis",label="X axis:", choices = columns_num, selected=columns_num[3])
                  ),
                  column(3,
                         textInput(inputId = "scatter_y_label", label = "Y axis label:", placeholder = "Custom label"),
                         selectInput("scatter_yaxis",label="Y axis:",choices = columns_num, selected=columns_num[4])
                  ),
                  column(3,
                         textInput(inputId = "scatter_z_label", label = "Color axis label:", placeholder = "Custom label"), 
                         selectInput(inputId = "scatter_zaxis", label = "Color axis:", choices = c("none", columns_num), selected=columns_num[5])
                         
                  ),
                  column(3,
                         selectInput(inputId = "scatter_color", label = "Color Type:", choices = c("lightgoldenrod1", "azure2"), selected = "lightgoldenrod1"),
                         checkboxInput(inputId = "scatter_round", label = "Round to Integer",value = F),
                         checkboxInput(inputId = "scatter_log10", label = "log10", value = F),
                         checkboxInput(inputId = "scatter_density", label = "density (bad performance)", value = F),
                         checkboxInput(inputId = "scatter_line", label = "line", value = T)
                  ),
                  fluidRow(
                    column(12,
                       actionButton(inputId = "scatter_plot", label="Plot",
                                    style = "color: #fff; background-color: #3c8dbc")
                )
              )
              )
            ),
            
            fluidRow(
              box(width=12, title="Output", id="Output_scatter",
                plotlyOutput("plot_scatter",height="100%")
              )
            )
    ),
    
    
    #Scatter_Category --------------------------------------------------------
    #scatter sub_page
    tabItem(tabName = "scatter_cat",
        fluidRow(
          box(width=12,title="Inputs",id="scatter_cat_inputs",
              column(3,
                     textInput(inputId = "scatter_cat_X_label", label = "X axis label:", placeholder = "Custom label"),
                     selectInput("scatter_cat_xaxis",label="X axis:", choices = columns_num, selected=columns_num[3])
              ),
              column(3,
                     textInput(inputId = "scatter_cat_y_label", label = "Y axis label:", placeholder = "Custom label"),
                     selectInput("scatter_cat_yaxis",label="Y axis:",choices = columns_num, selected=columns_num[4])
              ),
              column(3,
                     textInput(inputId = "scatter_cat_z_label", label = "Category label:", placeholder = "Custom label"), 
                     selectInput(inputId = "scatter_cat_zaxis", label = "Categories:", choices = c("none", columns), selected=columns[2])
                     
              ),
              column(3,
                     selectInput(inputId = "scatter_cat_color", label = "Color Type:", choices = c("lightgoldenrod1", "azure2"), selected = "lightgoldenrod1"),
                     checkboxInput(inputId = "scatter_cat_round", label = "Round to Integer",value = F),
                     checkboxInput(inputId = "scatter_cat_log10", label = "log10", value = F),
                     checkboxInput(inputId = "scatter_cat_density", label = "density (bad performance)", value = F),
                     checkboxInput(inputId = "scatter_cat_line", label = "line", value = T)
              ),
              fluidRow(
                column(12,
                   actionButton(inputId = "scatter_cat_plot", label="Plot",
                                style = "color: #fff; background-color: #3c8dbc")
            )
          )
          )
        ),
        
        fluidRow(
          box(width=12, title="Output", id="Output_scatter_cat",
              plotlyOutput("plot_scatter_cat",height="100%")
          )
        )
        ),




    
    #heatmap-----------------------------------------------------------------------
    #Main page for heatmap
    tabItem(tabName = "heatmap",
            fluidRow(
              box(width=12,title="Inputs",id="heatmap_inputs",
                  column(2,
                         selectInput("heat_select_row",label="Genes:",multiple=TRUE, choices = genes, selected=genes[c(2,3)]),
                    br(),
                         selectInput("heat_mode",label="Data transformation:",c("raw","log2","zscore"),selected="raw")
                  ),
                  column(2,
                         selectInput("heat_select_col", label="Columns:", multiple = T, choices = columns_num, selected = columns_num[c(1,2)]) #only numeric selectable
                  ),
                  column(2,
                         selectInput("heat_distrib",label="Data distribution",c("auto","one-sided","two-sided"), selected="auto"),
                    br(),
                         selectInput("heat_color",label="Color Type:",c("reds","viridis","plasma","inferno","magma", "blues", "heat", "cubehelix", "ylgnbu"),selected="BrBG")
                  ),
                  column(2,
                         selectInput("heat_clustering",label="Clustering",c("none", "row", "column", "both"), selected="both"), #both
                         br(),
                         selectInput("heat_clusterdist",label="Cluster Distance",c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"), selected="manhattan")
                  ),
                  column(2,
                         selectInput("heat_clustermethod",label="Cluster method",c("average", "ward.D", "ward.D2", "single", "complete", "mcquitty", "median", "centroid"), selected="average"),
                         br(),
                         textInput("heat_unitlabel",label="Unit label",value = "Enter unit...")
                  ),
                  column(1,
                         checkboxInput("heat_reverse",label="Reverse Coloring:",value = FALSE),
                         br(),
                         checkboxInput("heat_rowlabel",label="Row Label",value = T) #true
                  ),
                  column(1,
                         checkboxInput("heat_columnlabel",label="Column Label",value = T) #true
                  ),
                  fluidRow(
                    column(12,
                           actionButton(inputId = "heat_plot", label="Plot",
                                        style = "color: #fff; background-color: #3c8dbc")
                           )
                  )
              )
            ),
            
          
            #-u|--unitlabel       Label for unit displayed above color key.
            
            
            fluidRow(
              tabBox(width = NULL, height = NULL, title="Output", id="Output_heatmap",side="right", 
                     tabPanel("Plot", plotlyOutput("plot_heatmap", height = "auto", width = "auto")),
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
                  selectInput("gv_select",label="Genes:",multiple=TRUE, choices = genes, selected= genes[2])
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
                  ),
                  fluidRow(
                    column(12,
                       actionButton(inputId = "genview_plot", label="Plot",
                                    style = "color: #fff; background-color: #3c8dbc")
                )
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
  )
  )
)
    #-----------------------------------------------------------------------
  

# Server --------------------------------------------------------------------

server <- function(input, output, session) {  source("helpers.R") #for dev purpose
   #genview-------------------------------------------------------------
  #Section genview
  
  dataInput_genview <- reactive({
    genes <- table1[table1[[2]] %in% input$gv_select]
  })
  
  genview_table <- eventReactive(input$genview_plot, {
    dataInput_genview()
  })
  
  genview_plot <- eventReactive(input$genview_plot, {
    genes <- dataInput_genview()
    #only select numeric columns
    #values <- as.data.frame(genes[, sapply(genes, is.numeric), with = FALSE])
    #row.names(values) <- genes[[1]]

    #print(colnames(values))
    #values2<-values[,3:ncol(table1)]
    #width=1
    #Function Call for Genview Plots
    x <- dynamic_matrixsplit(genes, reps, input$gv_plottype, input$gv_facet,input$gv_color,input$gv_cols)#, width,input$gv_height)
    ggplotly(x,height=900, tooltip = "text")
  })
  
  output$table_genview<-renderDataTable({
    genview_table()
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
    genview_plot()
  })#, height=900)
  
  #heatmap-------------------------------------------------------------
  #Section heatmap
  
  #get selected data
  dataInput_heat <- reactive({
    #get selected rows (genes)
    genes<-table1[table1[[2]] %in% input$heat_select_row]
    #get selected columns
    cols <- c(colnames(table1)[1:2], input$heat_select_col)
    
    selectedData<- genes[, ..cols]
  })
  
  #change plot if button is pressed
  heatmap_plot <- eventReactive(input$heat_plot,{
    plot <- create_heatmaply(data = dataInput_heat(),
                             clustmethod = input$heat_clustermethod,
                             clustdist = input$heat_clusterdist,
                             clustering = input$heat_clustering,
                             color_vector = input$heat_color,
                             reverse_coloring = input$heat_reverse,
                             rowlabel = input$heat_rowlabel,
                             collabel = input$heat_columnlabel,
                             unitlabel = input$heat_unitlabel
                             )
    
    #print(plot)
    #plot_method "plotly" or row-dendrogram is rotated
    #plot <- heatmaply(plot, plot_method = "plotly") #%>%layout(margin = c(500,100,100,100))
    #print(attributes(plot))
    #plot
    #subplot(plot, plot)
  })
  
  #change table if button is pressed
  heatmap_table <- eventReactive(input$heat_plot, {
    dataInput_heat()
  })
  
  output$table_heatmap<-renderDataTable({
    heatmap_table()
  },options=list(scrollX=TRUE))
  
  output$plot_heatmap<-renderPlotly({
    heatmap_plot()
  })
  
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
  scatter_plot <- eventReactive(input$scatter_plot, {
    if(input$scatter_zaxis == "none"){
      selectedData <- table1[, c(colnames(table1)[1], input$scatter_xaxis, input$scatter_yaxis), with = F]
    }else{
      x <- c(colnames(table1)[1], input$scatter_xaxis, input$scatter_yaxis, input$scatter_zaxis)
      selectedData <- table1[, x, with = F]
    }
    
    #selectedData <- table1[, c(1,3,4,5)]
    #write.csv(selectedData,file="test.txt")
    #message(print(selectedData))
    ##print(str(selectedData))
    
    ##rows <- which(as.logical((selectedData[[2]]!=0) + (selectedData[[3]] != 0)))
    ##selectedData <- selectedData[rows]
    
    ##selectedData[,2:4] <- log10(selectedData[,2:4])
    ##is.na(selectedData) <- sapply(selectedData, is.infinite)
    ##selectedData[is.na(selectedData)] <- 0
    
    x <- create_scatterplot2(selectedData, input$scatter_round, input$scatter_log10, colors = input$scatter_color, x_label = input$scatter_X_label, y_label = input$scatter_y_label, z_label = input$scatter_z_label, density = input$scatter_density, line = input$scatter_line)
    ##x <- ggplot(selectedData, aes(selectedData[[2]], selectedData[[3]], color = selectedData[[4]])) + 
    ##  geom_point(alpha = 0.5) +
    ##  stat_density2d()
    ggplotly(x, tooltip = "all")
    #ggplot(table1,aes(x=input$scatter_xaxis,y=input$scatter_yaxis))# + geom_point() 
  })
  
  output$plot_scatter<-renderPlotly({
    scatter_plot()
  })#, height=550)
  
  # Section Scatter Category ------------------------------------------------
  scatter_cat_plot <- eventReactive(input$scatter_cat_plot, {
    if(input$scatter_cat_zaxis == "none"){
      selectedData <- table1[, c(colnames(table1)[1], input$scatter_cat_xaxis, input$scatter_cat_yaxis), with = F]
    }else{
      selectedData <- table1[, c(colnames(table1)[1], input$scatter_cat_xaxis, input$scatter_cat_yaxis, input$scatter_cat_zaxis), with = F]
    }
    
    x <- create_scatterplot2(selectedData, input$scatter_cat_round, input$scatter_cat_log10, colors = input$scatter_cat_color, x_label = input$scatter_cat_X_label, y_label = input$scatter_cat_y_label, z_label = input$scatter_cat_z_label, density = input$scatter_cat_density, line = input$scatter_cat_line, categorized = TRUE)
    ggplotly(x)
  })
  
  output$plot_scatter_cat<-renderPlotly({
    scatter_cat_plot()
  })#, height=550)
  
}




shinyApp(ui, server)