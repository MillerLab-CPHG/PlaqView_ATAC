###### PlaqView Master Code ######
##### Author: Wei Feng Ma, UVA.
##### wm5wt@virginia.edu

#### DATABASE NAMES AND COLOR SCHEMES ####
# below line is commented for shinyapp.io deployment temp
### set this once in terminal before deploying to shinyapps.io ###
# options(repos = BiocManager::repositories())

# enrichR functions 
# handcurate db names 
dbs <- c("KEGG_2019_Human",
         "WikiPathways_2019_Human",
         "GO_Biological_Process_2018",
         "ChEA_2016",
         "GWAS_Catalog_2019",
         "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
         "Gene_Perturbations_from_GEO_down",
         "Gene_Perturbations_from_GEO_up")
enrichRdb <- sort(dbs)

# color definitions
original_color_list <-
  {c("rosybrown2",
     "cadetblue1",
     "lemonchiffon3",
     "darkseagreen",
     "skyblue3",
     "thistle3",
     "cadetblue3",
     "darkseagreen1",
     "palevioletred3",
     "palevioletred1",
     "darkseagreen2",
     "rosybrown3",
     "thistle2",
     "lightsteelblue3",
     "salmon1",
     "palevioletred4",
     "lemonchiffon4",
     "cadetblue2"
  )}

color_function <- colorRampPalette(original_color_list)
# color_function <- colorRampPalette(metcolors)

manual_color_list <- color_function(40) # change this if clusters >40

#### LIBRARIES #### 
# library(BiocManager)
library(shiny)
library(shinythemes)
library(Seurat)
library(shinybusy) #install.packages("shinybusy")
library(enrichR) # install.packages("enrichR")
library(waiter)
library(DT)
library(readxl)
library(shinyWidgets)
library(shinyjs)
library(rDGIdb) # BiocManager::install("rDGIdb")
library(tidyverse)
library(ggpubr)
library(gtools)
library(ArchR)
library(hexbin) # this is required for archR
library(Signac)
library(parallel)
# library(reactlog)
# library(future)
#
# # tell shiny to log all reactivity
# reactlog_enable()
# 
# tell shiny to try to paralle compute
future::plan("multisession")
addArchRThreads(threads = 4) 


#### READ GOOGLE SHEET ####
googlesheets4::gs4_deauth() # this tells google sheet to read-only
df <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1hLyjPFA2ZRpBLHnTgUnmDz7kimMZWFbz_ZGTml3-hRA/edit#gid=0", sheet = "ATAC")


# df <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vTF5Gw4Dbshlh3wVB8UAMswUEiOn4NEzXaEp8x73NtbWY3n4oIrWEVNMIwNYyInJM7k70G1lUcr7x9g/pub?output=csv")

df$DOI <- paste("<a href=",  df$DOI,">", "Link", "</a>") # this converts to clickable format

# subset data rows that are marked 'deployed = Yes"
df <- filter(df, `Deployed` == "Yes")
df <- df %>% 
  select(Authors, Year, Journal, DOI, Species, Tissue, Notes, Population, Cell.Number, 'DataID', `Article.Title` ) 
df$`Article.Title` <- str_to_title(df$`Article.Title`) # autocaps

#### UI ####
# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Aesthetics ####
  theme = shinytheme("flatly"),
  add_busy_bar(color = "#ff9142", height = "100px"), # THIS IS THE BUSY BAR
  use_waiter(), 
  # waiter_show_on_load(html = spin_rotate()),
  useShinyjs(),
  # div(DT::dataTableOutput("table"), style = "font-size: 75%; width: 75%"), # DT font sinzes
  
  # Pages ####
  navbarPage("PlaqView-ATAC", id = "inTabset",
             
             #### UI: Data ####
             tabPanel("Select Dataset", 
                      mainPanel(width = 12,
                                fluidRow(
                                  column(width = 5,
                                         wellPanel(
                                           img(src = "logo.png", width = '100%'),
                                           h3("Instructions:"),
                                           tags$ol(
                                             tags$li("Start by clicking on desired dataset."),
                                             tags$li("Click the blue 'Load Dataset' button. (This button is disabled until you click on a dataset!)"), # change server code to toggle
                                             tags$li("The 'Start Exploring' button will appear when data is loaded."),
                                             tags$li("(Optional) come back to this page to load another dataset.")
                                             
                                           ),
                                           br(),
                                           fluidRow(
                                             column(width = 12,
                                                    # load data button
                                                    actionBttn(
                                                      inputId = "loaddatabutton",
                                                      label = "Load Dataset",
                                                      style = "unite",
                                                      color = "primary",
                                                      block = T),
                                                    
                                                    textOutput("loadeddatasetID"),
                                                    
                                                    br(),
                                                    # jump to page 1 button
                                                    hidden(
                                                      actionBttn(
                                                        inputId = "jumpto1",
                                                        label = "Start Exploring",
                                                        style = "unite",
                                                        color = "success",
                                                        block = T)
                                                    ),
                                                    
                                             ),
                                           
                                                    
                                             
                                           ),
                                           )
                                  ),
                                  column(width = 7,
                                         wellPanel(
                                           includeMarkdown("descriptionfiles/helptext_welcome.Rmd"),
                                           img(src = "abstract.png", width = '100%'),
                                         )
                                         ), # column
                                  
                                 
                                  
                                ) # fluid row 
                                ), # mainpanel
                      mainPanel(width = 12,
                                wellPanel(
                                  h3("Click to Select a Single-Cell Dataset"),
                                  br(),
                                  DT::dataTableOutput('availabledatasettable'),
                                          br(),

                                          ),
                              
                      )
             ), # tab panel
             
             #### UI: Genes   ----
             tabPanel(title = "Gene Lookup", value = "panel1",
                      mainPanel(width = 12, # 12/12 is full panel
                                fluidRow(## panel for gene input
                                  column(
                                    width = 5,
                                    wellPanel(                                      
                                      # must add up to 12 all columns
                                      textInput(
                                        "genes",
                                        width = '100%',
                                        h3("Query Gene Expression", h5("please follow HUGO conventions")),
                                        value = "NOX4",
                                        placeholder = "try: TREM2, CYBB"
                                      ),
                                      
                                      # choose the type of output graph 
                                      
                                      pickerInput(
                                        inputId = "selectlabelmethodforgenequery",
                                        label = "Select Labeling Method", 
                                        choices = list (
                                          "UMAP (Numbered) Clusters" = "Clusters",
                                          "Author Provided Annotation" = "Author_Provided"
                                        ), 
                                        selected = "UMAP Clusters",
                                        width = '95%' #neeed to fit this
                                      ),
                                      
                                      # 'go' button
                                      actionBttn(
                                        inputId = "runcode",
                                        label = "Query Gene Expression",
                                        style = "unite",
                                        color = "success",
                                        block = T)
                                      
                                    )
                                    
                                  ),
                                  
                                  ## panel for description
                                  column(
                                    width = 7,
                                    wellPanel(includeMarkdown("descriptionfiles/helptext_singlegenepage.Rmd"))
                                  )
                                ),
                                
                                
                                #spacer
                                br(),
                                
                                
                                ## lower panel for graphic outputs
                                wellPanel(width = 12,
                                          # textOutput("selecteddatasetID"),  
                                          fluidRow( # top split rows
                                            column(width = 6, align = "center", 
                                                   plotOutput("umaps", width = "auto", height = '500px'),
                                                   br(),
                                                   downloadButton("download.umap", "Download UMAP", width = '100%')
                                            ),
                                            
                                            column(width = 6, align="center", 
                                                   plotOutput("genequeryumap", width = "auto", height = '500px'),
                                                   br(),
                                                   downloadButton("download.genequeryumap", "Download Expression Plot", width = '100%')
                                                   
                                            ) # column 
                                            
                                          ), # fluidrow
                                          br(),
                                          fluidRow( # bottom whole for GO output
                                            column(width = 12, 
                                                   selectInput("selectedenrichRdb", label = h5("Top Significantly Enriched Pathways"), 
                                                               choices = enrichRdb, 
                                                               selected = "GO_Biological_Process_2018")
                                            ),
                                            column(width = 12, 
                                                   DT::dataTableOutput("enrichtable"),
                                                   helpText("You must restart query if you change database"),
                                                   downloadButton("downloadenrichRdata", "Download Complete Pathway Enrichment Data")
                                                   
                                            ),
                                          )# another fluidrow 
                                          
                                )# wellpanel
                                
                                
                      )# MAIN PANEL CLOSURE
             ), # TAB PANEL CLOSURE
             
             
             
             
             #### App UI Closure ####
  )# close navbarpage
  
)# close fluidpage



#### SERVER FUNCTIONS ####
# Define server logic required to draw a histogram
server <- function(input, output, session) {
  #### SER: Data ####
  output$availabledatasettable <-
    DT::renderDataTable(df, server = F, # server is for speed/loading
                        selection = list(mode = 'single'),
                        options=list(columnDefs = list(list(visible=FALSE, targets=c(10)))), # this hides the #8 col (datasetID)
                        escape = FALSE) # this escapes rendering html (link) literally and makes link clickable
  
  # start the page with load data disabled until dataset is clicked
  disable("loaddatabutton")
 
   observeEvent(input$loaddatabutton, {
    # create path for loading data
    path.to.archr.plaqviewobj <- file.path(paste("data/", df$DataID[input$availabledatasettable_rows_selected], sep=""))
  
    plaqviewobj <<- loadArchRProject(path = path.to.archr.plaqviewobj)
    
    # show which data is read
    loadeddatasetID <<- paste("Loaded Dataset: ", print(df$DataID[input$availabledatasettable_rows_selected]))
    output$loadeddatasetID <- renderText(loadeddatasetID)
    
    ## these are just for displaying current data name in other tabs##
    output$selecteddatasetID <- renderText({
      paste0("Current dataset: ", df$DataID[input$availabledatasettable_rows_selected])
    }) 
    output$selecteddatasetID2 <- renderText({
      paste0("Current dataset: ", df$DataID[input$availabledatasettable_rows_selected])
    }) 
    output$selecteddatasetID3 <- renderText({
      paste0("Current dataset: ", df$DataID[input$availabledatasettable_rows_selected])
    }) 
   

  })
   
   observeEvent(input$availabledatasettable_rows_selected,{
     enable("loaddatabutton")
   })

  # server code to show jumpto1
  observeEvent(input$loaddatabutton, {
    shinyjs::show(id = "jumpto1")  
  }) 
  
  # server code to jump to page 1
  observeEvent(input$jumpto1, {
    updateTabsetPanel(session = getDefaultReactiveDomain(), "inTabset",
                      selected = "panel1") # this is to switch to tab1
    
  })

  
  #### SER: Genes ####
  observeEvent(input$runcode,{ 
    
    #### NOMENCLATURE UPDATE
    if (df$Species[input$availabledatasettable_rows_selected] == "Human") {
      corrected <- str_to_upper(input$genes)
    } else{
      corrected <- str_to_title(input$genes)
    }
  
    parsed.genes <- str_split(input$genes, ", ")[[1]]
    
    updateTextInput(getDefaultReactiveDomain(),
                    "genes", value = corrected)
    #### UMAPS
    output$umaps <- 
      renderPlot(
        plotEmbedding(
          ArchRProj = plaqviewobj,
          colorBy = "cellColData",
          # UMAP
          name = input$selectlabelmethodforgenequery,
          # subheader
          embedding = "UMAP"
        )
      ) # render plot
       

    output$genequeryumap <- 
      renderPlot(
        plotEmbedding(
        ArchRProj = plaqviewobj, 
        colorBy = "GeneIntegrationMatrix", 
        name = parsed.genes, 
        embedding = "UMAP",
        imputeWeights = getImputeWeights(plaqviewobj)
      )
      )# render plot
        
      
    #### DOWNLOAD UMAPS ####
    #### download label umap 
    output$download.umap<- downloadHandler(
      filename = function() {
        paste(df$DataID[input$availabledatasettable_rows_selected], "_", 
              "UMAP.pdf", sep = "")
      },
      content = function(file) {
        pdf(file, onefile = F)
        temp <- plotEmbedding(
          ArchRProj = plaqviewobj,
          colorBy = "cellColData",
          # UMAP
          name = input$selectlabelmethodforgenequery,
          # subheader
          embedding = "UMAP"
        )
        
        plot(temp)
        dev.off()
      }
    )# close downloadhandler
    
    #### download gene umap 
    output$download.genequeryumap<- downloadHandler(
      filename = function() {
        paste(df$DataID[input$availabledatasettable_rows_selected], "_", input$genes,
              "UMAP.pdf", sep = "")
      },
      content = function(file) {
        pdf(file, onefile = F)
        
        temp <- plotEmbedding(
          ArchRProj = plaqviewobj, 
          colorBy = "GeneIntegrationMatrix", 
          name = parsed.genes, 
          embedding = "UMAP",
          imputeWeights = getImputeWeights(plaqviewobj)
        )
        
        plot(temp)
        dev.off()
      }
      
    )# close downloadhandler
    
    
    #### GSEA ####
    enriched <- enrichr(genes = parsed.genes, 
                        database = enrichRdb) # this queries all of them
    cleanedenrichedtable <- select(enriched[[input$selectedenrichRdb]], -Old.Adjusted.P.value, -Old.P.value,)
    cleanedenrichedtable <- top_n(cleanedenrichedtable, 100) # top 100 will be rendered
    
    #select columns to display
    cleanedenrichedtable <- cleanedenrichedtable %>% select(Term, Overlap, Adjusted.P.value, Combined.Score, Genes)
    
    # force as.numeric to remove a bug in DT pkg
    #cleanedenrichedtable$Adjusted.P.value <- as.numeric(cleanedenrichedtable$Adjusted.P.value)
    #cleanedenrichedtable$Adjusted.P.value <- as.numeric(cleanedenrichedtable$Combined.Score)
    
    output$enrichtable <- DT::renderDataTable(cleanedenrichedtable, server = F)
    
    # Downloadable csv of selected dataset
    output$downloadenrichRdata <- downloadHandler(
      filename = function() {
        paste(parsed.genes, "_pathwayenrichment.csv", sep = "")
      },
      content = function(file) {
        write.csv(enriched[[input$selectedenrichRdb]], file, row.names = FALSE)
      } 
    )# close downloadhandler
    

    
  })# closes observe event
  

  
  

  
  
} # ends server function

# Run the application 
shinyApp(ui = ui, server = server)

