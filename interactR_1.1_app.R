library(shiny)
library(dplyr)
library(DT)
library(magrittr)

# custom function
`%notin%` = function(x,y) !(x %in% y)

# load the full table data set
df = readRDS("data/InteractR_interactome_full.Rds")

# interpro domain annotations for human proteome
interpro = readRDS("data/Interpro_for_interactR.Rds")

# data types to exclude
to_exc =
  c("Two-hybrid","Proximity Label-MS","IP-MS","Affinity Capture-Luminescence","Affinity Capture-Western","Biochemical Activity",
    "Co-crystal Structure","Far Western","FRET",
    "Reconstituted Complex","Co-localization","Protein-peptide",
    "PCA","Affinity Capture-RNA","Protein-RNA",
    "Co-purification","Co-fractionation")

# Define UI ----
ui <- fluidPage(
  
  
  
  titlePanel(tags$h1("InteractR")),
  tags$div(
    tags$h3("Exploring protein interactors"),
    tags$h6("v1.1 Dominic D.G.Owens (dominic.owens at utoronto.ca)")
  ),
  
  sidebarLayout(
    sidebarPanel(
      #h3("Explore protein interactors"),
      
      # search bar
      textInput(inputId = "insearch",
                label = p("HGNC gene symbol"), 
                placeholder ="e.g. TP53",
                value = ""),
      
      #p("View all interactions with metadata or non-redundant interactors"),
      # full of small table
      radioButtons("display", h3("Options"),
                   choices = list("Full table" = 2, "Non-redundant interactors (with domains)" = 1),selected = 2),
      
      # exclude data types
      checkboxGroupInput("data_to_exclude",
                         label = "Select data types to exclude",
                         choices = to_exc),
      
      #p("View or download result"),
      # download
      downloadButton("download", "Download filtered data")),
    
    mainPanel(
      h3("Results"),
      DTOutput(outputId = "table"),
      h4("Data from the following sources:"),
      p("STRING v12.0 physical links network, experimental evidence only (not text mining) https://string-db.org/"),
      p("Human cell map (BFDR threshold=0.1) https://humancellmap.org/"),
      p("BioGRID v4.4.222 physical interactions https://thebiogrid.org"),
      p("BioPlex (IP-MS) https://bioplex.hms.harvard.edu/"),
      p("OpenCell (IP-MS) https://opencell.czbiohub.org/download"),
      p("Human Reference interactome (Y2H) http://www.interactome-atlas.org/"),
      p("HIPPIE v2.3 http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/download.php)"),
      p("Interpro human domains annotation https://www.ebi.ac.uk/interpro/"),
      p("updated June 28 2023")
    )
  )
  
  
)


# Define server logic ----
server <- function(input, output) {
  
  # Filter data based on search term
  filtered_df <- reactive({
    # require inputs before displaying to avoid errors
    req(input$insearch, input$display) 
    if (input$insearch == "") {
      p("Search results will be displayed here")
      # return the full table if radioButtons set to full table (2)
    } else if (input$display == 2) {
      #Sys.sleep(1)
      
      # convert gene names to upper case
      this_input = toupper(input$insearch)
      
      # exclude the data types we don't want
      ## note we also need to exclude the "unknown" data types
      data_types_to_exc = c(input$data_to_exclude,"unknown")
      df %<>%
        filter(experiment_type %notin% data_types_to_exc)
      
      #df = df[(df[, 1] %notin% data_types_to_exc)]
      
      # filter actually on this gene
      df[(df[, 1] == this_input)|(df[, 3] == this_input), ]
      
      # return the non-redundant interactors table if radioButtons set (1)
    } else if (input$display == 1){
      #Sys.sleep(1)
      
      # convert gene names to uppercase
      this_input = toupper(input$insearch) 
      
      # exclude the data types we don't want
      ## note we also need to exclude the "unknown" data types
      data_types_to_exc = c(input$data_to_exclude,"unknown")
      df %<>%
        filter(experiment_type %notin% data_types_to_exc)
      #df = df[(df[, 1] %in% data_types_to_exc)]
      
      
      a=df[(df[, 1] == this_input)|(df[, 3] == this_input), ]
      ## clean up to just have the interactors
      a = unique(c(a$A_Gene,a$B_Gene))
      # remove itself from the list
      a=sort(a[a!=this_input])
      
      # make a tibble and give nice colname
      a = as_tibble(a)
      colnames(a) = "Gene"
      
      # join to the interpro domains file
      a %<>% left_join(interpro,by="Gene")
      a
    }
  })
  
  # Display filtered data in table
  output$table <- renderDT({
    filtered_df()
  })
  
  
  # Download filtered data as TSV file
  output$download <- downloadHandler(
    filename = function() {
      paste("filtered_data_", toupper(input$insearch), ".tsv", sep = "")
    },
    content = function(file) {
      write.table(filtered_df(), file, sep = "\t", row.names = FALSE)
    }
  )
  
  
  
  
}

# Run the app ----
shinyApp(ui = ui, server = server)
