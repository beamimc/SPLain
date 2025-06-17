# launch app.R
# This script is used to launch a Shiny app with the required data files.

library(here)
library(shiny)

# Ensure the working directory is set to the project root
setwd(here::here()) 

# Load the required data files
se <- readRDS(here::here("data", "glinos_saturn_dtu.rds"))
exons <- readRDS(here::here("data", "glinos_exons.rds"))

# Load the app.R file from the app directory
app_dir <- file.path(here::here(), "app")

# somewhere mention all dependencies and how to install e.g.
# pkgs <- renv::dependencies("app/app.R")
# find.package(pkgs$Package)

# but also "clusterProfiler"?

# "BSgenome.Hsapiens.UCSC.hg38" --> also mention hg38 somewhere as the default genome
source(file.path(app_dir, "app.R"))

# Run the Shiny app with the provided data
shiny::runApp(app(
  se = se,
  exons = exons,
  app_dir = app_dir
))
