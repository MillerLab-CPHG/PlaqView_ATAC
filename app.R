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
library(reactable)

# library(reactlog)
# library(future)
#
# # tell shiny to log all reactivity
# reactlog_enable()
# 
# tell shiny to try to paralle compute
# future::plan("multisession")
addArchRThreads(threads = 1) 


#### READ GOOGLE SHEET ####
googlesheets4::gs4_deauth() # this tells google sheet to read-only
df <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1hLyjPFA2ZRpBLHnTgUnmDz7kimMZWFbz_ZGTml3-hRA/edit#gid=0", sheet = "ATAC")


# df <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vTF5Gw4Dbshlh3wVB8UAMswUEiOn4NEzXaEp8x73NtbWY3n4oIrWEVNMIwNYyInJM7k70G1lUcr7x9g/pub?output=csv")

df$DOI <- paste("<a href=",  df$DOI,">", "Link", "</a>") # this converts to clickable format

# subset data rows that are marked 'deployed = Yes"
df <- filter(df, `Deployed` == "Yes")
df <- df %>% 
  select('DataID', Year, Journal, DOI, Species, Tissue, Notes, Population, Cells = Cell.Number, `Article.Title` ) 
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
                                             tags$li("Select a dataset from the BLUE drop-down Menu."),
                                             tags$li("Click the blue 'Load Dataset' button."), # change server code to toggle
                                             tags$li("The 'Start Exploring' button will appear when data is loaded."),
                                             tags$li("(Optional) come back to this page to load another dataset.")
                                             
                                           ),
                                           br(),
                                           fluidRow(
                                             column(width = 12,
                                                    # load data button
                                                    h4("Select a Dataset"),
                                                    pickerInput(
                                                      inputId = "dataselector",
                                                      #label = "Select a Dataset", 
                                                      choices = df$DataID,
                                                      choicesOpt = list(
                                                        subtext = paste(df$Species, ": ",
                                                                        df$Cells,
                                                                        " Cells",
                                                                        sep = "")),
                                                      options = list(
                                                        style = "btn-primary")
                                                    ),
                                                    
                                                    actionBttn(
                                                      inputId = "loaddatabutton",
                                                      label = "Load Dataset",
                                                      style = "unite",
                                                      color = "primary",
                                                      block = T),
                                                    
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
                                                    br(),
                                                    helpText(textOutput("loadeddatasetID")) ,
                                                    
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
                                
                                h4("Details of Single- Cell Dataset and IDs"),
                                actionButton(inputId = "refreshtable", "Fetch Latest Dataset Details"),
                                reactableOutput("availabledatasettable"),
                                br(),
                                inlineCSS(list("table" = "font-size: 12px")),
                                
                                
                                
                      )
             ),
             
             #### UI: Genes   ----
             tabPanel(title = "Gene Score Predictor", value = "panel1",
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
             
             
             
             
             #### UI: Track   ----
             tabPanel(title = "Track Browser", value = "panel1",
                      mainPanel(width = 12, # 12/12 is full panel
                                fluidRow(## panel for gene input
                                  column(
                                    width = 5,
                                    wellPanel(                                      
                                      # must add up to 12 all columns
                                      textInput(
                                        "genesfortrack",
                                        width = '100%',
                                        h3("Query Chromatin Tracks", h5("please follow HUGO conventions")),
                                        value = "NOX4",
                                        placeholder = "try: TREM2, CYBB"
                                      ),
                                      
                                      # choose the type of output graph 
                                      
                                      pickerInput(
                                        inputId = "selectlabelmethodfortrackquery",
                                        label = "Select Labeling Method", 
                                        choices = list (
                                          "UMAP (Numbered) Clusters" = "Clusters",
                                          "Author Provided Annotation" = "Author_Provided"
                                        ), 
                                        selected = "UMAP Clusters",
                                        width = '95%' #neeed to fit this
                                      ),
                                      materialSwitch(
                                        inputId = "peak2geneswitch",
                                        label = "Peak2Gene Overlay", 
                                        value = FALSE,
                                        status = "primary"
                                      ),
                                      # 'go' button
                                      actionBttn(
                                        inputId = "runtrack",
                                        label = "Query Genome Track",
                                        style = "unite",
                                        color = "success",
                                        block = T)
                                      
                                    )
                                    
                                  ),
                                  
                                  ## panel for description
                                  column(
                                    width = 7,
                                    wellPanel(includeMarkdown("descriptionfiles/helptext_trackpage.Rmd"))
                                  )
                                ),
                                
                                
                                #spacer
                                br(),
                                
                                
                                ## lower panel for graphic outputs
                                wellPanel(width = 12,
                                          fluidRow(
                                            column(width = 12, 
                                                   plotOutput("genometrack", width = "auto", height = '1000px'),
                                                   
                                                   )
                                          )
                                          
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
    renderReactable({
      reactable(df, compact = T, searchable = T, defaultPageSize = 20,
                defaultColDef = colDef(
                  header = function(value) gsub(".", " ", value, fixed = TRUE),
                  # ell = function(value) format(value, nsmall = 1),
                  align = "center",
                  minWidth  = 50,
                  headerStyle = list(background = "#f7f7f8")
                ),
                columns = list(
                  'Article.Title' = colDef(minWidth = 200,
                  ),
                  # DOI = colDef(html = TRUE, cell = function(value, index) {
                  #   # this is raw html
                  #   sprintf('<a href="%s" target="_blank">LINK</a>', df$DOI[index], value)
                  #   }),
                  Cells = colDef(format = colFormat(separators = TRUE),
                                 footer = paste0("Total Cells: ", sum(as.numeric(df$Cells)))
                  ),
                  
                  Journal = colDef(html = TRUE, cell = function(value, index) {
                    sprintf('<a href="%s" target="_blank">%s</a>', df$DOI[index], value)
                  }),       
                  Year = colDef(show = T),
                  DOI = colDef(show = F),
                  
                  Tissue = colDef(show = F),
                  Notes = colDef(show = F),
                  Species = colDef(show = T),
                  DataID = colDef(
                    minWidth = 100,
                    name = "Data ID",
                    html = TRUE,
                    cell = function(value, index){
                      species <- df$Species[index]
                      year <- df$Year[index]
                      tissue <- df$Tissue[index]
                      Notes <- df$Notes[index]
                      DOI <- df$DOI[index]
                      div(
                        div(style = list(fontWeight = 600), value),
                        #div(style = list(fontSize = 12), species),
                        #div(style = list(fontSize = 12), year),
                        div(style = list(fontSize = 12), tissue),
                        div(style = list(fontSize = 12), Notes)
                      )# div
                    }
                  )
                )
      )
      
      
    })
  
  # refresh button 
  
  # refresh button
  observeEvent(input$refreshtable, {
    session$reload()
    session$reload()
    session$reload()
    session$reload()
    session$reload()
    session$reload()
    
    
  })
  
  observeEvent(input$loaddatabutton, {
    # create path for loading data
    path <- file.path(paste("data/", input$dataselector, "/", sep=""))
    plaqviewobj <<- loadArchRProject(path)

    # show which data is read
    loadeddatasetID <<- paste("Dataset Loaded Sucessfully: ", print(input$dataselector))
    output$loadeddatasetID <- renderText(loadeddatasetID)
    
    ## these are just for displaying current data name in other tabs##
    output$selecteddatasetID <- renderText({
      paste0("Current dataset: ", input$dataselector)
    }) 
    output$selecteddatasetID2 <- renderText({
      paste0("Current dataset: ", input$dataselector)
    }) 
    output$selecteddatasetID3 <- renderText({
      paste0("Current dataset: ", input$dataselector)
    }) 
    
    
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
    if (df$Species[df$DataID == input$dataselector] == "Human") {
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
        paste(input$dataselector, "_", 
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
        paste(input$dataselector, "_", input$genes,
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
        paste(parsed.genes, input$dataselector, "_pathwayenrichment.csv", sep = "")
      },
      content = function(file) {
        write.csv(enriched[[input$selectedenrichRdb]], file, row.names = FALSE)
      } 
    )# close downloadhandler
    

    
  })# closes observe event
  

  
  

  
  
  #### SER: Track ####
  observeEvent(input$runtrack,{ 
    
    #### NOMENCLATURE UPDATE
    if (df$Species[df$DataID == input$dataselector] == "Human") {
      corrected <- str_to_upper(input$genesfortrack)
    } else{
      corrected <- str_to_title(input$genesfortrack)
    }
    
    parsed.genes <- str_split(input$genesfortrack, ", ")[[1]]
    
    updateTextInput(getDefaultReactiveDomain(),
                    "genesfortrack", value = corrected)

    output$genometrack <- 
      renderPlot({
        if(input$peak2geneswitch == TRUE){
          
          
          
          temp <- plotBrowserTrack(
            ArchRProj = plaqviewobj, 
            groupBy = input$selectlabelmethodfortrackquery, 
            geneSymbol = input$genesfortrack, 
            upstream = 50000,
            downstream = 50000, 
            loops = getPeak2GeneLinks(plaqviewobj)
          )
          

          
        }else{
          

          temp <- plotBrowserTrack(
            ArchRProj = plaqviewobj, 
            groupBy = input$selectlabelmethodfortrackquery, 
            geneSymbol = input$genesfortrack, 
            upstream = 50000,
            downstream = 50000 
          )
        }
        
        
        grid::grid.newpage()
        grid::grid.draw(temp[[input$genesfortrack]])
          
          }
      ) # render plot
    
    #### DOWNLOAD TRACK ####
    #### download label umap 
    output$download.umap<- downloadHandler(
      filename = function() {
        paste(input$dataselector, "_", 
              "Genome_Browser_Tracks.pdf", sep = "")
      },
      content = function(file) {
        pdf()
        
  
      }
    )# close downloadhandler
    

    
    
    
  })# closes observe event
  
  
  
  
  
  
  
  
} # ends server function

# Run the application 
shinyApp(ui = ui, server = server)

