isoformAnalysisUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidPage(

      # Global DT table font-size styling and fixed smaller height for first row
      tags$style(HTML("
        .dataTables_wrapper table {
          font-size: 20px;
        }
        .small-row {
          height: 250px;  /* Adjust height as needed */
           auto;
        }
      ")),

      # ─── First Row (smaller height) ─────────────────────────────
      fluidRow(
        class = "border-bottom p-3",   # <- added small-row class here
        # DTU table (left)
        column(
          width = 6,
          tags$h4("Significant DTU transcripts"),
          DTOutput(ns("dtu_table"))
        ),
        # Gene description + GO (right)
        column(
          width = 6,
          tags$div(
            # strong("Gene Description:"),
            textOutput(ns("gene_description"))
          ),
          tabsetPanel(
            id   = ns("go_tabs"),
            type = "tabs",
            tabPanel("GO Terms", DTOutput(ns("go_table")))
          )
        )
      ),

      # ─── Second Row (larger height by default) ──────────────────
      fluidRow(
        class = "border-bottom p-3",

        # Isoform Plot card (6 cols)
        column(
          width = 6,
          tags$div(class = "card mb-3",
                   tags$div(class = "card-header bg-primary text-white", "Isoform Structures"),
                   tags$div(class = "card-body",
                            plotlyOutput(ns("isoform_plot"), width = "100%")
                   )
          )
        ),

        # Mean Differences card (3 cols)
        column(
          width = 3,
          tags$div(class = "card mb-3",
                   tags$div(class = "card-header bg-primary text-white", 
                            textOutput(ns("barplot_title"))),
                   tags$div(class = "card-body",
                            plotlyOutput(ns("barplot"), width = "100%")
                   )
          )
        ),

        # Transcript Proportions card (3 cols)
        column(
          width = 3,
          tags$div(class = "card mb-3",
                   tags$div(class = "card-header bg-primary text-white", 
                            "Transcript Proportions by Condition"),
                   tags$div(class = "card-body",
                            plotlyOutput(ns("lineplot"), width = "100%")
                   )
          )
        )
      )

    ) # fluidPage
  ) # tagList
}



# Module server functions
isoformAnalysisServer <- function(id, se, exons, dtu_df, sig_res, selected_conditions) {
  moduleServer(id, function(input, output, session) {


    output$barplot_title <- renderText({
      cd1 <- selected_conditions()$cd1
      cd2 <- selected_conditions()$cd2
      conds <- selected_conditions()
      paste0("Mean Differences between Conditions (", cd2, " - ", cd1, ")")
    })
    
    
    # no updateSelectizeInput needed
    
    selected_gene <- reactive({
      sel <- input$dtu_table_rows_selected
      if (!is.null(sel) && length(sel) == 1) {
        dtu_df()[["Symbol"]][sel]
      } else {
        # default to row 1 when nothing is selected
        dtu_df()[["Symbol"]][1]
      }
    })
    
    
    # Reactive values
    mean_diffs_DTU <- reactive({
      req(selected_gene())
      calc_mean_diff_DTU(se, selected_gene(), sig_res() )
    })
    
    prop <- reactive({
      req(selected_gene())
      calc_prop(se, selected_gene(), sig_res())
    })
    
    pvals <- reactive({
      req(selected_gene())
      get_pvals(se, selected_gene(), sig_res())
    })
    
    # plot_path <- reactive({
    #   req(selected_gene())
    #   temp_file <- "www/temp.png"
    #   plot_gene_txs(selected_gene(), temp_file, mean_diffs_DTU(), pvals())
    #   temp_file
    # })
    
    # Outputs
    output$gene_description <- renderText({
      req(selected_gene())
      gene_info <- get_description(selected_gene())
      gene_info <- paste0("Gene description: ", gene_info)
      if (length(gene_info) > 0) gene_info else "Select a gene to see description."
    })
    
    output$dtu_table <- renderDT({
      datatable(dtu_df(), selection = "single", options = list(pageLength = 3))
    })
    
    output$go_table <- renderDT({
      req(selected_gene())
      go_data <- get_GO(selected_gene())
      datatable(go_data, options = list(pageLength = 2))
    })
    
    output$go_plot <- renderPlot({
      req(selected_gene())
      go_plot(c(selected_gene()))
      
    })
    
    # output$gene_plot <- renderImage({
    #   req(plot_path())
    #   list(src = plot_path(), contentType = "image/png", width = "100%", height = "auto")
    # }, deleteFile = FALSE)
    # 
    output$isoform_plot <- renderPlotly({
      req(selected_gene())
      plot_isoforms_wiggle(exons, selected_gene(), sig_res())
    })
    
    # Example: before calling plotting functions
    txp_colors <- reactive({
      transcripts <- unique(c(colnames(prop()), names(mean_diffs_DTU())))
      get_transcript_colors(transcripts)
    })
    
    
    output$barplot <- renderPlotly({
      req(mean_diffs_DTU(), pvals())
      barplot_meandifs(mean_diffs_DTU(), pvals(), txp_colors())
      
      
    })
    
    output$lineplot <- renderPlotly({
      req(prop())
      cd1 <- selected_conditions()$cd1
      cd2 <- selected_conditions()$cd2
      line_plot_txp_comparison(se, prop(), txp_colors(),cd1, cd2 )
    })
  })
}
