library(shiny)
library(bslib)
library(shinycssloaders)
library(DT)
library(plotly)
library(reactable)
library(dplyr)
library(ggplot2)
library(scales)
library(purrr)
library(stringr)
library(grid)
library(patchwork)

library(SummarizedExperiment)
library(GenomicFeatures)
library(GenomicRanges)
library(org.Hs.eg.db)
library(GO.db)
library(biomaRt)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)

library(clusterProfiler)
library(enrichplot)
library(plyranges)
library(reshape2)
library(viridis)
library(RColorBrewer)
library(here)

# # Load Local Development Package -----------------------------------------------
# setwd("/work/users/b/e/beacm/wiggleplotr")
# devtools::load_all()

# Set Working Directory and Source Scripts -------------------------------------
setwd(here::here())  
wd <- getwd() 

app <- function(se, exons, app_dir = ".") {
  shiny::addResourcePath(
    prefix = "assets",
    directoryPath = file.path(app_dir, "www")
  )
  
  #load functions in app/R/
  helpers_env <- new.env()
  r_scripts <- list.files(file.path(app_dir, "R"), pattern = "\\.R$", full.names = TRUE)
  invisible(lapply(r_scripts, function(f) source(f, local = helpers_env)))
  list2env(as.list(helpers_env), envir = environment())
  
  #load data
  data <- list(se = se, exons = exons)
  se <- data$se 
  
  condition_choices_df <- parse_saturnDTU_conditions(se)
  data$condition_choices <-  condition_choices_df
  
  
  source(file.path(app_dir, "ui.R"), local = TRUE)
  ui <- build_ui(condition_choices)
  server_func <- source(file.path(app_dir, "server.R"), local = TRUE)$value
  
  #launch App
  shinyApp(
    ui = ui,
    server = function(input, output, session) {
      server_func(input, output, session, data)
    }
  )
}
