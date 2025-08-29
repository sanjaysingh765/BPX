# --- Load Required Libraries ---
library(shiny)
library(bslib)
library(ggplot2)
library(shinyWidgets)
library(shinycssloaders)
library(DT)
library(readxl)
library(purrr)
library(dplyr)
library(tidyr)
library(stringr)

# --- Shiny UI Definition ---
ui <- page_fluid(
  theme = bs_theme(
    version = 5,
    bootswatch = "cerulean",
    base_font = font_google("Roboto"),
    heading_font = font_google("Montserrat"),
    bg = "#f8f9fa",
    fg = "#212529",
    primary = "#0073e6"
  ),
  
  # App Title
  titlePanel("BioRad Advanced qPCR Analyzer"),
  
  layout_sidebar(
    sidebar = sidebar(
      title = "Controls",
      h4("Instructions"),
      helpText("1. Upload the required Excel files.", 
               "2. Customize and analyze plots."),
      hr(),
      
      # File inputs
      fileInput("amplification_file", "Upload Amplification Data (.xlsx)",
                accept = ".xlsx", buttonLabel = "Browse...", placeholder = "No file selected"),
      
      fileInput("plate_layout_file", "Upload Plate Layout (.xlsx)",
                accept = ".xlsx", buttonLabel = "Browse...", placeholder = "No file selected"),
      
      hr(),
      uiOutput("fluorophore_selector_ui"),
      uiOutput("sample_selector_ui"),
      
      hr(),
      h4("Plot Customization"),
      materialSwitch("show_average", "Show average of replicates", status = "primary"),
      materialSwitch("show_cq", "Show Cq lines", status = "info"),
      radioGroupButtons("plot_display", "Display Mode", 
                        choices = c("All Together", "Faceted by Sample"), justified = TRUE),
      
      sliderInput("font_size", "Plot Font Size", 8, 22, 14),
      sliderInput("line_width", "Line Width", 0.2, 2, 0.8, step = 0.1),
      
      conditionalPanel(
        condition = "input.show_cq == true && input.tabs == 'Amplification Plot'",
        hr(),
        sliderInput("cq_font_size", "Cq Label Font Size", 1, 8, 3.5, step = 0.5)
      ),
      
      hr(),
      downloadBttn("download_data", "Download Plotted Data", style = "fill", color = "primary")
    ),
    
    # --- Main Content ---
    mainPanel(
      tabsetPanel(
        id = "tabs",
        type = "pills",
        
        # Welcome Page
        tabPanel(
          "Welcome",
          card(
            full_screen = FALSE,
            card_header(
              tags$img(src = "https://biopathogenix.com/wp-content/uploads/2021/07/BioPathogenix-Horizontal-1.svg",
                       style = "max-width: 60%; height: auto;")
            ),
            card_body(
              h3("Welcome to the BioRad Advanced qPCR Analyzer"),
              p("This application provides a powerful and intuitive interface for analyzing 
                BioRad qPCR amplification data. You can:"),
              tags$ul(
                tags$li("Upload Amplification Data and Plate Layout (.xlsx)"),
                tags$li("Visualize amplification plots with customizable options"),
                tags$li("Overlay replicate averages and Cq threshold lines"),
                tags$li("Perform precision analysis of Cq values"),
                tags$li("Export plots and processed data for reporting")
              ),
              p("Navigate through the tabs to begin analysis.")
            )
          )
        ),
        
        # Amplification Plot
        tabPanel(
          "Amplification Plot",
          card(
            full_screen = TRUE,
            card_header("Amplification Curves"),
            card_body(
              withSpinner(plotOutput("amplification_plot", height = "700px"), 
                          type = 6, color = "#0073e6")
            ),
            card_footer(
              downloadBttn("download_amp_plot", "Download Plot", style = "bordered", color = "success")
            )
          )
        ),
        
        # Precision Analysis
        tabPanel(
          "Precision Analysis",
          card(
            full_screen = TRUE,
            card_header("Precision of Cq Values"),
            card_body(
              withSpinner(plotOutput("precision_plot", height = "500px"), type = 6),
              hr(),
              withSpinner(DTOutput("precision_table"), type = 6)
            ),
            card_footer(
              downloadBttn("download_precision_plot", "Download Plot", style = "bordered", color = "success"),
              downloadBttn("download_precision_table", "Download Table", style = "bordered", color = "primary")
            )
          )
        )
      ),
      
      # Footer Credit
      hr(),
      div(
        style = "text-align: center; color: #555; font-size: 14px; padding: 10px;",
        HTML("&copy; 2025 <b>Biopathogenix</b>. All rights reserved.")
      )
    )
  )
)

# --- Shiny Server Logic ---
server <- function(input, output, session) {
  
  # --- DATA LOADING AND JOINING ---
  
  plate_layout_with_cq <- reactive({
    req(input$plate_layout_file)
    df <- read_excel(input$plate_layout_file$datapath, sheet = 1)
    validate(need(all(c("Well", "Sample", "Cq") %in% names(df)), 
                  "Plate layout must contain 'Well', 'Sample', and 'Cq' columns."))
    df %>% 
      mutate(
        Well = toupper(str_trim(Well)),
        # FIX: Changed [A-H] to [A-Z] to support 384-well plates and beyond.
        Row = str_extract(Well, "^[A-Z]"), 
        Col = str_extract(Well, "[0-9]+$"),
        Well_formatted = paste0(Row, sprintf("%02d", as.numeric(Col))),
        Cq = as.numeric(Cq)
      ) %>%
      select(Well_formatted, Sample, Cq)
  })
  
  amplification_data <- reactive({
    req(input$amplification_file)
    path <- input$amplification_file$datapath
    data_sheets <- excel_sheets(path) %>% setdiff("Run Information")
    map_dfr(data_sheets, function(sheet) {
      read_excel(path, sheet = sheet) %>%
        pivot_longer(cols = -Cycle, names_to = "Well", values_to = "RFU") %>%
        mutate(Fluorophore = sheet)
    })
  })
  
  all_data_joined <- reactive({
    req(amplification_data(), plate_layout_with_cq())
    amplification_data() %>%
      mutate(
        Well = toupper(str_trim(Well)),
        # FIX: Changed [A-H] to [A-Z] here as well for consistency.
        Row = str_extract(Well, "^[A-Z]"), 
        Col = str_extract(Well, "[0-9]+$"),
        Well_formatted = paste0(Row, sprintf("%02d", as.numeric(Col)))
      ) %>%
      left_join(plate_layout_with_cq(), by = "Well_formatted") %>%
      mutate(Sample = ifelse(is.na(Sample), "Unassigned", Sample))
  })
  
  # --- DYNAMIC UI CONTROLS ---
  
  output$fluorophore_selector_ui <- renderUI({
    req(all_data_joined())
    fluorophores <- unique(all_data_joined()$Fluorophore)
    selectInput("selected_fluorophore", "Select Fluorophore", choices = fluorophores)
  })
  
  output$sample_selector_ui <- renderUI({
    req(all_data_joined())
    samples <- unique(all_data_joined()$Sample)
    checkboxGroupInput("selected_samples", "Select Samples", choices = samples, selected = samples)
  })
  
  # --- FILTERING AND AVERAGING ---
  
  plot_data_filtered <- reactive({
    req(all_data_joined(), input$selected_fluorophore, input$selected_samples)
    all_data_joined() %>%
      filter(Fluorophore == input$selected_fluorophore, Sample %in% input$selected_samples)
  })
  
  plot_data_averaged <- reactive({
    req(plot_data_filtered())
    plot_data_filtered() %>%
      group_by(Sample, Cycle) %>%
      summarise(
        RFU_mean = mean(RFU, na.rm = TRUE), RFU_sd = sd(RFU, na.rm = TRUE),
        Cq_mean = mean(Cq, na.rm = TRUE), .groups = "drop"
      )
  })
  
  # --- PRECISION ANALYSIS ---
  
  precision_data <- reactive({
    req(plot_data_filtered())
    plot_data_filtered() %>%
      filter(!is.na(Cq)) %>%
      distinct(Well, Sample, Cq) %>%
      group_by(Sample) %>%
      summarise(
        Mean_Cq = mean(Cq, na.rm = TRUE),
        SD_Cq = sd(Cq, na.rm = TRUE),
        CV_Cq = (SD_Cq / Mean_Cq) * 100,
        .groups = "drop"
      )
  })
  
  precision_plot_object <- reactive({
    req(plot_data_filtered())
    
    df <- plot_data_filtered() %>% filter(!is.na(Cq)) %>% distinct(Well, Sample, Cq)
    
    ggplot(df, aes(x = Sample, y = Cq, color = Sample)) +
      geom_jitter(width = 0.2, size = 3, alpha = 0.8) +
      stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "black") +
      theme_minimal(base_size = input$font_size) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none") +
      labs(title = "Precision of Cq Values by Sample",
           x = "Sample",
           y = "Cq Value")
  })
  
  output$precision_plot <- renderPlot({
    precision_plot_object()
  })
  
  output$precision_table <- renderDT({
    datatable(precision_data(), options = list(pageLength = 10), rownames = FALSE) %>%
      formatRound(columns = c('Mean_Cq', 'SD_Cq', 'CV_Cq'), digits = 2)
  })
  
  # --- AMPLIFICATION PLOT ---
  
  amplification_plot_object <- reactive({
    req(plot_data_filtered())
    y_max <- max(plot_data_filtered()$RFU, na.rm = TRUE)
    p <- ggplot()
    
    if (input$show_average) {
      df_avg <- plot_data_averaged()
      p <- p + geom_line(data = df_avg, aes(x = Cycle, y = RFU_mean, color = Sample), linewidth = input$line_width) +
        geom_ribbon(data = df_avg, aes(x = Cycle, ymin = RFU_mean - RFU_sd, ymax = RFU_mean + RFU_sd, fill = Sample), alpha = 0.2, color = NA)
      if(input$show_cq) {
        cq_label_data <- df_avg %>% filter(!is.na(Cq_mean)) %>% distinct(Sample, Cq_mean)
        p <- p + geom_vline(data = cq_label_data, aes(xintercept = Cq_mean, color = Sample), linetype = "dashed", linewidth = input$line_width * 0.8) +
          geom_text(data = cq_label_data, aes(x = Cq_mean, y = y_max, label = round(Cq_mean, 2), color = Sample), size = input$cq_font_size, show.legend = FALSE, angle = 90, vjust = 1.5, hjust = 1)
      }
    } else {
      df_ind <- plot_data_filtered()
      p <- p + geom_line(data = df_ind, aes(x = Cycle, y = RFU, group = Well, color = Sample), alpha = 0.8, linewidth = input$line_width)
      if(input$show_cq) {
        cq_lines_data <- df_ind %>% filter(!is.na(Cq)) %>% distinct(Well, Sample, Cq)
        cq_label_data <- df_ind %>% filter(!is.na(Cq)) %>% group_by(Sample) %>% summarise(Cq_mean = mean(Cq, na.rm = TRUE), .groups = "drop")
        p <- p + geom_vline(data = cq_lines_data, aes(xintercept = Cq, color = Sample), linetype = "dashed", linewidth = input$line_width * 0.8) +
          geom_text(data = cq_label_data, aes(x = Cq_mean, y = y_max, label = round(Cq_mean, 2), color = Sample), size = input$cq_font_size, show.legend = FALSE, angle = 90, vjust = 1.5, hjust = 1)
      }
    }
    
    p <- p + theme_minimal(base_size = input$font_size) +
      labs(title = paste("Amplification Plot for", input$selected_fluorophore), 
           x = "Cycle", y = "Relative Fluorescence Units (RFU)", 
           color = "Sample", fill = "Sample") +
      theme(legend.text = element_text(size = rel(0.9)))
    
    if (input$plot_display == "Faceted by Sample") { p <- p + facet_wrap(~ Sample) }
    
    return(p)
  })
  
  output$amplification_plot <- renderPlot({
    amplification_plot_object()
  })
  
  # --- DOWNLOAD HANDLERS ---
  
  output$download_data <- downloadHandler(
    filename = function() {
      paste0("plotted_data-", input$selected_fluorophore, "-", Sys.Date(), ".csv")
    },
    content = function(file) {
      df <- if (input$show_average) plot_data_averaged() else plot_data_filtered()
      write.csv(df, file, row.names = FALSE)
    }
  )
  
  output$download_amp_plot <- downloadHandler(
    filename = function() { paste0("amplification_plot-", input$selected_fluorophore, "-", Sys.Date(), ".png") },
    content = function(file) {
      ggsave(file, plot = amplification_plot_object(), device = "png", width = 12, height = 8, dpi = 300)
    }
  )
  
  output$download_precision_plot <- downloadHandler(
    filename = function() { paste0("precision_plot-", input$selected_fluorophore, "-", Sys.Date(), ".png") },
    content = function(file) {
      ggsave(file, plot = precision_plot_object(), device = "png", width = 10, height = 7, dpi = 300)
    }
  )
  
  output$download_precision_table <- downloadHandler(
    filename = function() { paste0("precision_data-", input$selected_fluorophore, "-", Sys.Date(), ".csv") },
    content = function(file) {
      write.csv(precision_data(), file, row.names = FALSE)
    }
  )
}

# --- Run the Application ---
shinyApp(ui = ui, server = server)
