#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(dplyr)
library(ggplot2)

# Load datasets
rna_data <- read.csv("~/Desktop/HPA/rna_tissue_consensus.tsv", sep = "\t", stringsAsFactors = FALSE)
scRNA_data <- read.csv("~/Desktop/HPA/rna_single_cell_type.tsv", sep = "\t", stringsAsFactors = FALSE)
ihc_data <- read.csv("~/Desktop/HPA/normal_ihc_data.tsv", sep = "\t", stringsAsFactors = FALSE)

# Define UI
ui <- fluidPage(
  titlePanel("GeneExpresso: Gene Expression Analysis using Human Protein Atlas(HPA)"),
  
  sidebarLayout(
    sidebarPanel(
      textInput("gene", "Enter Gene Name:", value = "CR2"),
      numericInput("fold_consensus", "Fold Change Threshold (Consensus Data):", value = 3, min = 1),
      numericInput("fold_single_cell", "Fold Change Threshold (Single-Cell Data):", value = 2, min = 1),
      actionButton("analyze", "Analyze")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Tissue Expression (Consensus Data)", plotOutput("consensusPlot")),
        tabPanel("Cell Type Expression (Single-Cell Data)", plotOutput("singleCellPlot")),
        tabPanel("Tissue Selectivity Results", tableOutput("consensusTable")),
        tabPanel("Cell Type Selectivity Results", tableOutput("singleCellTable")),
        tabPanel("IHC Reliability", tableOutput("ihcTable"))
      )
    )
  )
)

# Define Server
server <- function(input, output) {
  observeEvent(input$analyze, {
    gene_name <- input$gene
    fold_change_threshold_consensus <- input$fold_consensus
    fold_change_threshold_single_cell <- input$fold_single_cell
    
    # Filter and visualize consensus data
    output$consensusPlot <- renderPlot({
      gene_data <- rna_data %>% filter(Gene.name == gene_name)
      if (nrow(gene_data) > 0) {
        ggplot(gene_data, aes(x = Tissue, y = nTPM, fill = Tissue)) +
          geom_bar(stat = "identity") +
          labs(title = paste("Tissue Expression for Gene:", gene_name),
               x = "Tissue",
               y = "nTPM") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
      } else {
        plot.new()
        title("No data found for this gene in consensus RNA data.")
      }
    })
    
    # Filter and visualize single-cell data
    output$singleCellPlot <- renderPlot({
      cell_data <- scRNA_data %>% filter(Gene.name == gene_name)
      if (nrow(cell_data) > 0) {
        ggplot(cell_data, aes(x = Cell.type, y = nTPM, fill = Cell.type)) +
          geom_bar(stat = "identity") +
          labs(title = paste("Cell Type Expression for Gene:", gene_name),
               x = "Cell Type",
               y = "nTPM") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
      } else {
        plot.new()
        title("No data found for this gene in single-cell RNA data.")
      }
    })
    
    # Calculate tissue selectivity (consensus data)
    output$consensusTable <- renderTable({
      gene_data <- rna_data %>% filter(Gene.name == gene_name)
      if (nrow(gene_data) > 0) {
        overall_mean <- mean(gene_data$nTPM, na.rm = TRUE)
        gene_data %>%
          mutate(Fold_Change = nTPM / overall_mean) %>%
          filter(Fold_Change >= fold_change_threshold_consensus)
      } else {
        data.frame(Message = "No data found for this gene in consensus RNA data.")
      }
    })
    
    # Calculate cell type selectivity (single-cell data)
    output$singleCellTable <- renderTable({
      cell_data <- scRNA_data %>% filter(Gene.name == gene_name)
      if (nrow(cell_data) > 0) {
        overall_mean <- mean(cell_data$nTPM, na.rm = TRUE)
        cell_data %>%
          mutate(Fold_Change = nTPM / overall_mean) %>%
          filter(Fold_Change >= fold_change_threshold_single_cell)
      } else {
        data.frame(Message = "No data found for this gene in single-cell RNA data.")
      }
    })
    
    # Fetch IHC reliability
    output$ihcTable <- renderTable({
      ihc_data %>%
        filter(Gene.name == gene_name) %>%
        select(Gene, Tissue, Reliability)
    })
  })
}

# Run App
shinyApp(ui = ui, server = server)
