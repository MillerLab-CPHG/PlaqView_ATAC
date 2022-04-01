temp <- plotBrowserTrack(
  ArchRProj = plaqviewobj, 
  groupBy = "Author_Provided", 
  geneSymbol = "NOX4", 
  upstream = 50000,
  downstream = 50000 
)



grid::grid.newpage()
grid::grid.draw(temp[["NOX4"]])
