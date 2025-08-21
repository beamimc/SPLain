exonLevelUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      class = "border-bottom p-3",
      
      column(width = 6,
        tags$h4("Significant DTU transcripts"),
        DTOutput(ns("dtu_table"))
      ),
      column(width = 6,
             tags$div(class = "card mb-3",
                      tags$div(class = "card-header bg-primary text-white", 
                               "Isoform Structures with selected exons"),
                      tags$div(class = "card-body",
                               plotlyOutput(ns("exon_level_plot"), width = "100%"))
        )
      )
    ),
    # Tabs for downstream windows and another plot
    tabsetPanel(
      id = ns("exon_tabs"),  # optional id
      tabPanel(
        "Window Plots",
        fluidRow(
          column(
            width = 6,
            tags$div(class = "card mb-3",
                     tags$div(class = "card-header bg-primary text-white", 
                              "Up and downstream sequences"),
                     tags$div(class = "card-body",
                              plotOutput(ns("window_summary_plot_down"), width = "100%"))
            )
          ),
          column(
            width = 6,
            tags$div(class = "card mb-3",
                     tags$div(class = "card-header bg-primary text-white", 
                              "Up and downstream sequences"),
                     tags$div(class = "card-body",
                              plotOutput(ns("window_summary_plot_nonreg"), width = "100%"))
            )
          )
        )
      ),
      tabPanel(
        "Single Plot",
        fluidRow(
          class = "border-bottom p-3",
          column(
            width = 12,
            tags$div(class = "card mb-3",
                     tags$div(class = "card-header bg-primary text-white", 
                              "Alternative View"),
                     tags$div(class = "card-body",
                              plotOutput(ns("plot_window_comparison"), width = "100%"))
            )
          )
        )
      )
    )
  )
}


exonLevelServer <- function(id, exons, dtu_df, x_flat, sig_res, ref_assembly) {
  moduleServer(id, function(input, output, session) {
    selected_gene <- reactive({
  req(filtered_dtu())
  df <- filtered_dtu()
  sel <- input$dtu_table_rows_selected
  idx <- if (!is.null(sel) && length(sel) == 1) sel else 1
  df[["Symbol"]][idx]
})

    

    output$dtu_table <- DT::renderDT({
      df <- filtered_dtu()
      DT::datatable(df, selection = "single", options = list(pageLength = 5))
    })

    # Reactive: Downregulated exons and windows
    downstream_data <- reactive({
      req(x_flat())
      downreg_exons <- detect_downreg_exonsv2(x_flat())
      validate(need(length(downreg_exons) > 0, "No downregulated exons found"))
      list(
        exons = downreg_exons,
        downstream = get_downstream_from_GRanges(downreg_exons, ref_assembly = ref_assembly),
        upstream = get_upstream_from_GRanges(downreg_exons, ref_assembly = ref_assembly)
      )
    })

    # Reactive: Nonregulated exons and windows
    nonreg_data <- reactive({
      req(x_flat())
      downreg_exons <- detect_downreg_exonsv2(x_flat())
      nonreg_exons <- get_nonreg_exons(x_flat(), downreg_exons)
      validate(need(length(nonreg_exons) > 0, "No nonregulated exons found"))
      list(
        exons = nonreg_exons,
        downstream = get_downstream_from_GRanges(nonreg_exons, ref_assembly = ref_assembly),
        upstream = get_upstream_from_GRanges(nonreg_exons, ref_assembly = ref_assembly)
      )
    })

    filtered_dtu <- reactive({
      req(downstream_data(), dtu_df())
      filter_txps <- names(downstream_data()$exons)
      df <- dtu_df()
      validate(need("Transcript" %in% names(df), "Column 'Transcript' not found in dtu_df()"))
      out <- df[df[["Transcript"]] %in% filter_txps, , drop = FALSE]
      validate(need(nrow(out) > 0, "No transcripts with downregulated exons"))
      out
    })

    

    output$exon_level_plot <- renderPlotly({
      req(selected_gene(), downstream_data())
      plot_downreg_exons(exons, selected_gene(), sig_res(), downstream_data()$exons)
    })

    output$window_summary_plot_down <- renderPlot({
      req(downstream_data())
      plot_updownstream_windows(
        downstream_data()$upstream,
        downstream_data()$downstream,
        exon_label = "downreg"
      )
    })

    output$window_summary_plot_nonreg <- renderPlot({
      req(nonreg_data())
      plot_updownstream_windows(
        nonreg_data()$upstream,
        nonreg_data()$downstream,
        exon_label = "non-reg"
      )
    })

    output$plot_window_comparison <- renderPlot({
      req(nonreg_data(), 
          downstream_data())
      plot_window_comparison(
        downstream_data()$upstream,
        nonreg_data()$upstream
      )
    })
  })
}
