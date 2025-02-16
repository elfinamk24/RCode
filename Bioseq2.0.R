# BioSeqAnalyzer: An Enhanced R Bioinformatics App with Colorful UI

# Load required packages
library(shiny)
library(Biostrings)
library(DT)
library(DECIPHER)

# Define UI
ui <- fluidPage(
  tags$head(
    tags$style(HTML("\n      body { background-color: #f7f7f7; font-family: Arial, sans-serif; }\n      .title-panel { background-color: #0073e6; color: white; padding: 15px; text-align: center; border-radius: 10px; }\n      .sidebar { background-color: #e6f2ff; padding: 15px; border-radius: 10px; }\n      .main-panel { background-color: #ffffff; padding: 20px; border-radius: 10px; box-shadow: 2px 2px 10px rgba(0, 0, 0, 0.1); }\n      .btn-primary { background-color: #0073e6; border: none; }\n      .btn-danger { background-color: #e60000; border: none; }\n    "))
  ),
  div(class = "title-panel", titlePanel("BioSeqAnalyzer: Sequence Analysis Tool by El <3")),
  sidebarLayout(
    sidebarPanel(
      class = "sidebar",
      fileInput("file", "Upload FASTA File", accept = ".fasta"),
      textInput("sequence", "Or Enter Sequence"),
      actionButton("analyze", "Analyze", class = "btn-primary"),
      actionButton("translate", "Translate to Protein", class = "btn-primary"),
      actionButton("align", "Perform Sequence Alignment", class = "btn-danger"),
      textInput("motif", "Enter Motif for Search")
    ),
    mainPanel(
      class = "main-panel",
      h3("Sequence Analysis"),
      DTOutput("seq_table"),
      plotOutput("gc_plot"),
      verbatimTextOutput("motif_search"),
      verbatimTextOutput("protein_translation"),
      verbatimTextOutput("alignment_result")
    )
  )
)

# Define server logic
server <- function(input, output) {
  
  analyze_seq <- eventReactive(input$analyze, {
    if (!is.null(input$file)) {
      seq_data <- readDNAStringSet(input$file$datapath)
    } else {
      seq_data <- DNAStringSet(input$sequence)
    }
    
    seq_length <- width(seq_data)
    gc_content <- letterFrequency(seq_data, c("G", "C"), as.prob = TRUE)
    gc_percent <- rowSums(gc_content) * 100
    
    data.frame(Sequence = names(seq_data), Length = seq_length, GC_Content = gc_percent)
  })
  
  output$seq_table <- renderDT({
    analyze_seq()
  })
  
  output$gc_plot <- renderPlot({
    df <- analyze_seq()
    barplot(df$GC_Content, names.arg = df$Sequence, col = "#0073e6", main = "GC Content (%)")
  })
  
  output$motif_search <- renderPrint({
    req(input$motif)
    seqs <- analyze_seq()
    motif_hits <- vmatchPattern(input$motif, DNAStringSet(seqs$Sequence))
    paste("Motif Found at Positions:", motif_hits)
  })
  
  output$protein_translation <- renderPrint({
    seqs <- analyze_seq()
    translated <- translate(DNAStringSet(seqs$Sequence))
    paste("Protein Sequence:", as.character(translated))
  })
  
  output$alignment_result <- renderPrint({
    if (!is.null(input$file)) {
      seq_data <- readDNAStringSet(input$file$datapath)
      aln <- AlignSeqs(seq_data)
      paste("Multiple Sequence Alignment:\n", aln)
    } else {
      "No sequences available for alignment."
    }
  })
}

# Run the application
shinyApp(ui = ui, server = server)
