build_ui <- function(condition_choices) {
    # at the very top of ui.R
  message("📍 Working directory is: ", normalizePath("."))
  message("📂 Files here: ", paste(list.files(), collapse = ", "))
  message("📂 www/: ", paste(list.files("www"), collapse = ", "))

  page_sidebar(
    title = tags$span(
    # 1. Your custom icon, sized & margin-adjusted
    tags$img(
      src    = "assets/logo_white.png",     # must live in www/icon_splain.png
      height = "50px",           # adjust to taste
      style  = "vertical-align: middle; margin-right: 8px;"
    ),
    # 2. Your app title
    "SPLain",
    style = "font-size:2rem; font-weight:bold; vertical-align: middle;"
  ),
      
    theme   = bs_theme(version = 5, bootswatch = "cosmo",    primary    = "#3d65a1"  ), 
       # Global CSS to increase overall font size
    tags$style(HTML("
      body, .sidebar, .card, .nav-tabs, .form-label, h1, h2, h3, h4, h5, h6 {
        font-size: 22px !important;
      }
      .card-header {
        font-size: 22px !important;
      }

    ")),


    sidebar = sidebarPanel(
  width = 12,

  # pick the DTU result column
  selectInput(
    inputId = "dtu_column",
    label   = "DTU comparison column",
    choices = "(loading...)"  # filled in server
  ),
  textOutput("current_column"),
  br(),

  fluidRow(
    column(
      width = 6,
      radioButtons(
        inputId  = "cd1",
        label    = "Reference",
        choices  = c("Pick a column first" = "loading_placeholder")
      )
    ),
    column(
      width = 6,
      radioButtons(
        inputId  = "cd2",
        label    = "Contrast",
        choices  = c("Pick a reference first" = "loading_placeholder")
      )
    )
  ),

  actionButton("apply_pair", "Apply comparison"),
  textOutput("current_comparison"),
      

  br(), br(),

  sliderInput("fdr_threshold", "FDR threshold:",
              min = 0, max = 1, value = 0.05, step = 0.01),

  actionButton("apply_fdr", "Apply FDR Filter"),
  textOutput("current_fdr"),
      

  br(), br(),

  radioButtons(
    inputId = "exon_filter",
    label = "Select structural change:",
    choices = c("Downregulated exon", "Upregulated exon", "UTR change"),
    selected = "Downregulated exon"
  )
    ),
    
    tabsetPanel(
      id = "main_tabs",
      tabPanel("DTU Exploration", value = "dtu", isoformAnalysisUI("isoform")),
      tabPanel("Exon-Level", value = "exon", exonLevelUI("exon")),
      tabPanel("Summary", value = "summary", summaryStatsUI("summary"))
    )
  )
}