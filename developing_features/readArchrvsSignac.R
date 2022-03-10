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
# library(RColorBrewer)
library(rDGIdb) # BiocManager::install("rDGIdb")
library(tidyverse)
# library(rsconnect)
library(monocle3)
library(ggpubr)
library(gtools)
library(CIPR)
library(ArchR)
library(Signac)
library(parallel)

# library(reactlog)
# library(future)
#
# # tell shiny to log all reactivity
# reactlog_enable()
# 
# # tell shiny to try to paralle compute
# future::plan("multisession")



# time code
ptm <- proc.time()
addArchRThreads(threads = 1) # doesnt seem like archr supports multithreading here

#### ArchR ####
proj <- loadArchRProject(path = "data/Turner_2022_ATAC/")


plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")


# proj <- addImputeWeights(proj)

plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneIntegrationMatrix", 
  name = "MYH11", 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(proj)
)


plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneIntegrationMatrix", 
  name = "MYH11", 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(proj)
)




p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "Clusters", 
  geneSymbol = "MYH11", 
  upstream = 50000,
  downstream = 50000
)


grid::grid.newpage()
grid::grid.draw(p$MYH11)

# output time code
proc.time() - ptm

#### Convert to Signac ####
