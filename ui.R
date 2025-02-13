# Define FELIXapp UI 
ui <- fluidPage(
  # Set UI theme 
  theme = bs_theme(version = 5, bootswatch = "minty"),
  # Sidebar layout with a input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      img(src='oatman2.jpg', height = "100%", width = "100%", align = "top"),
      
      # Input1: Upload r/s pathotype file ----
      fileInput(inputId = "df",
                label = "Choose r/s file",
                accept = c(".csv")),
      
      # Input2: Upload raw pathotype file ----
      fileInput(inputId = "df_raw",
                label = "Choose raw scores file",
                accept = c(".csv")),
      
      # Input3: Selector for choosing output type ----
      selectInput(inputId = "format",
                  label = "Choose output format:",
                  choices = c("Gene postulations (R/S)", "Gene postulations (Raw scores)", "Full sorted dataset (R/S)", "Full sorted dataset (Raw scores)")),
      
      #Input4: Choose the cultivar
      selectInput(inputId = "cult",
                  label = "Choose a cultivar",
                  choices = "",
                  selected = ""),
      
      # Input5: Numeric entry for number of obs to view ----
      numericInput(inputId = "obs",
                   label = "Number of cultivars to view:",
                   value = 25),
      # Download button
      downloadButton("downloadData", "Download")
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: HTML table with requested number of observations ----
      tableOutput("view"),
      
    )
  )
)