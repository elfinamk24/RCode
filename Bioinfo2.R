# Step 1: Load Libraries
library(shiny)
library(bio3d)
library(shinyWidgets)
library(DT)
library(shinythemes)
library(r3dmol)

# Step 2: Define UI
ui <- fluidPage(
  theme = shinytheme("cosmo"),
  titlePanel("Elfina Beats Pymol Huahaha"),
  sidebarLayout(
    sidebarPanel(
      textInput("pdb_id", "Enter PDB ID:", value = "1CRN"),
      actionButton("fetch", "Fetch Structure"),
      helpText("Example PDB IDs: 1CRN, 1CTF, 1IGT")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Summary", verbatimTextOutput("summary")),
        tabPanel("Atom Table", DTOutput("atom_table")),
        tabPanel("3D Structure", r3dmolOutput("pdb_view"))
      )
    )
  )
)

# Step 3: Define Server
server <- function(input, output) {
  pdb_data <- reactiveVal()
  
  observeEvent(input$fetch, {
    req(input$pdb_id)
    pdb <- tryCatch(read.pdb(input$pdb_id), error = function(e) NULL)
    if (!is.null(pdb)) {
      pdb_data(pdb)
    } else {
      showNotification("Invalid PDB ID or fetch error. Try another ID.", type = "error")
    }
  })
  
  output$summary <- renderPrint({
    req(pdb_data())
    cat("PDB ID:", input$pdb_id, "\n")
    cat("Protein Name:", pdb_data()$header$title, "\n")
    cat("Organism:", pdb_data()$header$organism, "\n")
  })
  
  output$atom_table <- renderDT({
    req(pdb_data())
    datatable(pdb_data()$atom, options = list(pageLength = 10))
  })
  
  output$pdb_view <- renderR3dmol({
    req(pdb_data())
    
    # Write the PDB file in the correct format
    write.pdb(pdb_data(), file = "temp.pdb")
    
    # Read the PDB file back in
    pdb_trimmed <- read.pdb("temp.pdb")
    
    # Convert PDB structure to a character string
    pdb_text <- paste0(readLines(pdb_trimmed$file), collapse = "\n")
    
    r3dmol() %>%
      m_add_model(pdb_text, format = "pdb") %>%
      m_set_style(style = m_style_cartoon()) %>%
      m_zoom_to()
  })
  
}

# Step 4: Run the App
shinyApp(ui = ui, server = server)