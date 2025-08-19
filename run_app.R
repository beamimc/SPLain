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

# "BSgenome.Hsapiens.UCSC.hg38" --> also mention hg38 somewhere as the default genome
source(file.path(app_dir, "app.R"))

# Run the Shiny app with the provided data
shiny::runApp(app(
  ref_assembly = "hg38",
  se = se,
  exons = exons,
  app_dir = app_dir
))
