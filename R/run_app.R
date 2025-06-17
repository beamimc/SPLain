# launch app.R
# This script is used to launch a Shiny app with the required data files.

# library(here)
# library(shiny)

# # Ensure the working directory is set to the project root
# setwd(here::here())

# # Load the required data files
# see <- readRDS(here::here("data", "glinos_saturn_dtu.rds"))
# exonss <- readRDS(here::here("data", "glinos_exons.rds"))

# Load the app.R file from the app directory
# app_dir <-file.path(here::here(), "app")
# source(file.path(app_dir, "app.R"))
# # Run the Shiny app with the provided data
# shiny::runApp(app(se= see,
#                   exons=exonss,
#                   app_dir = app_dir))

#' @export
run_SPLain_app <- function(se = NULL, exons = NULL, ...) {
  se <- se %||% readRDS(app_sys("data", "glinos_saturn_dtu.rds"))
  exons <- exons %||% readRDS(app_sys("data", "glinos_exons.rds"))
  app <- SPLain_app(se = se, exons = exons, app_dir = app_sys("app"))
  shiny::runApp(app, ...)
}
