################################################################################
################################  PACKAGES #####################################
################################################################################

library(shiny)
library(bslib)


################################################################################
###############################  Functions #####################################
################################################################################

rawscores_sorted <- function(df_raw, df, cult) {
  
  
  ##### List all cultivars in df and rank them by similarity to cultivar "cult"
  
  # Add a matches column to the first column spot in df 
  df <- cbind(Matches = rep(0, nrow(df)), df)
  
  #calculate the number of matches between cult and each cultivar
  for (i in 1:nrow(df)) {
    df$Matches[i] <- 0
    for (j in 4:ncol(df)) {
      if (is.na(df[i, j]) | is.na(df[df$Cultivar == cult, j])) {
        df$Matches[i] <- df$Matches[i] + 0
      }
      else if (df[i, j] == df[df$Cultivar == cult, j]) {
        df$Matches[i] <- df$Matches[i] + 1
      }
    }
  }
  
  # Add the df$Matches column to the first column spot in df_faw
  df_raw <- cbind(df$Matches, df_raw)
  colnames(df_raw)[1] <- "Matches"
  
  ##### List of cultivars with any that contain gene that can't be in cult removed
  
  df_trim <- df
  tocut <- c()
  
  for (c in 4:ncol(df_trim)){
    
    for (r in 1:nrow(df_trim)){
      
      # If the variety of interest is susceptible to a pathotype and the differential is not susceptible  #
      # to that pathotype, remove that differential/gene from the list of candidate genes                 #
      if(isTRUE(df_trim[df_trim$Cultivar == cult, c] == "s" && df_trim[r,c] == "r") == TRUE){
        tocut <- append(tocut, df_trim[r,2])
      }}}
  
  ##### Remove varieties containing impossible genes, stored in tocut, from df_trim
  
  # Trim tocut to only unique values and remove those from df_trim
  tocut <- unique(tocut)
  df_trim <- df_trim[!(df_trim$Cultivar %in% tocut),]
  
  #### Sort df, df_raw, df_trim and df_trim_raw by matches and prepare for output
  
  #sort df_trim by matches
  df_trim <- df_trim[order(df_trim$Matches, decreasing = TRUE),]
  
  #creat df_trim_raw with only the df_raw rows that are in df_trim
  df_trim_raw <- df_raw[df_raw$Cultivar %in% df_trim$Cultivar,]
  
  #sort df_trim_raw by matches
  df_trim_raw <- df_trim_raw[order(df_trim_raw$Matches, decreasing = TRUE),]
  
  #sort df by matches
  df <- df[order(df$Matches, decreasing = TRUE),]
  
  #sort df_rAW by matches
  df_raw <- df_raw[order(df_raw$Matches, decreasing = TRUE),]
  
}
rsscores_sorted <- function(df_raw, df, cult) {
  
  
  ##### List all cultivars in df and rank them by similarity to cultivar "cult"
  
  # Add a matches column to the first column spot in df 
  df <- cbind(Matches = rep(0, nrow(df)), df)
  
  #calculate the number of matches between cult and each cultivar
  for (i in 1:nrow(df)) {
    df$Matches[i] <- 0
    for (j in 4:ncol(df)) {
      if (is.na(df[i, j]) | is.na(df[df$Cultivar == cult, j])) {
        df$Matches[i] <- df$Matches[i] + 0
      }
      else if (df[i, j] == df[df$Cultivar == cult, j]) {
        df$Matches[i] <- df$Matches[i] + 1
      }
    }
  }
  
  # Add the df$Matches column to the first column spot in df_faw
  df_raw <- cbind(df$Matches, df_raw)
  colnames(df_raw)[1] <- "Matches"
  
  ##### List of cultivars with any that contain gene that can't be in cult removed
  
  df_trim <- df
  tocut <- c()
  
  for (c in 4:ncol(df_trim)){
    
    for (r in 1:nrow(df_trim)){
      
      # If the variety of interest is susceptible to a pathotype and the differential is not susceptible  #
      # to that pathotype, remove that differential/gene from the list of candidate genes                 #
      if(isTRUE(df_trim[df_trim$Cultivar == cult, c] == "s" && df_trim[r,c] == "r") == TRUE){
        tocut <- append(tocut, df_trim[r,2])
      }}}
  
  ##### Remove varieties containing impossible genes, stored in tocut, from df_trim
  
  # Trim tocut to only unique values and remove those from df_trim
  tocut <- unique(tocut)
  df_trim <- df_trim[!(df_trim$Cultivar %in% tocut),]
  
  #### Sort df, df_raw, df_trim and df_trim_raw by matches and prepare for output
  
  #sort df_trim by matches
  df_trim <- df_trim[order(df_trim$Matches, decreasing = TRUE),]
  
  #creat df_trim_raw with only the df_raw rows that are in df_trim
  df_trim_raw <- df_raw[df_raw$Cultivar %in% df_trim$Cultivar,]
  
  #sort df_trim_raw by matches
  df_trim_raw <- df_trim_raw[order(df_trim_raw$Matches, decreasing = TRUE),]
  
  #sort df by matches
  df <- df[order(df$Matches, decreasing = TRUE),]
  
}
rawscores_trimmed <- function(df_raw, df, cult) {
  
  
  ##### List all cultivars in df and rank them by similarity to cultivar "cult"
  
  # Add a matches column to the first column spot in df 
  df <- cbind(Matches = rep(0, nrow(df)), df)
  
  #calculate the number of matches between cult and each cultivar
  for (i in 1:nrow(df)) {
    df$Matches[i] <- 0
    for (j in 4:ncol(df)) {
      if (is.na(df[i, j]) | is.na(df[df$Cultivar == cult, j])) {
        df$Matches[i] <- df$Matches[i] + 0
      }
      else if (df[i, j] == df[df$Cultivar == cult, j]) {
        df$Matches[i] <- df$Matches[i] + 1
      }
    }
  }
  
  # Add the df$Matches column to the first column spot in df_faw
  df_raw <- cbind(df$Matches, df_raw)
  colnames(df_raw)[1] <- "Matches"
  
  ##### List of cultivars with any that contain gene that can't be in cult removed
  
  df_trim <- df
  tocut <- c()
  
  for (c in 4:ncol(df_trim)){
    
    for (r in 1:nrow(df_trim)){
      
      # If the variety of interest is susceptible to a pathotype and the differential is not susceptible  #
      # to that pathotype, remove that differential/gene from the list of candidate genes                 #
      if(isTRUE(df_trim[df_trim$Cultivar == cult, c] == "s" && df_trim[r,c] == "r") == TRUE){
        tocut <- append(tocut, df_trim[r,2])
      }}}
  
  ##### Remove varieties containing impossible genes, stored in tocut, from df_trim
  
  # Trim tocut to only unique values and remove those from df_trim
  tocut <- unique(tocut)
  df_trim <- df_trim[!(df_trim$Cultivar %in% tocut),]
  
  #### Sort df, df_raw, df_trim and df_trim_raw by matches and prepare for output
  
  #creat df_trim_raw with only the df_raw rows that are in df_trim
  df_trim_raw <- df_raw[df_raw$Cultivar %in% df_trim$Cultivar,]
  
  #sort df_trim by matches
  df_trim_raw <- df_trim_raw[order(df_trim_raw$Matches, decreasing = TRUE),]
  
}
rsscores_trimmed <- function(df_raw, df, cult) {
  
  
  ##### List all cultivars in df and rank them by similarity to cultivar "cult"
  
  # Add a matches column to the first column spot in df 
  df <- cbind(Matches = rep(0, nrow(df)), df)
  
  #calculate the number of matches between cult and each cultivar
  for (i in 1:nrow(df)) {
    df$Matches[i] <- 0
    for (j in 4:ncol(df)) {
      if (is.na(df[i, j]) | is.na(df[df$Cultivar == cult, j])) {
        df$Matches[i] <- df$Matches[i] + 0
      }
      else if (df[i, j] == df[df$Cultivar == cult, j]) {
        df$Matches[i] <- df$Matches[i] + 1
      }
    }
  }
  
  # Add the df$Matches column to the first column spot in df_faw
  df_raw <- cbind(df$Matches, df_raw)
  colnames(df_raw)[1] <- "Matches"
  
  ##### List of cultivars with any that contain gene that can't be in cult removed
  
  df_trim <- df
  tocut <- c()
  
  for (c in 4:ncol(df_trim)){
    
    for (r in 1:nrow(df_trim)){
      
      # If the variety of interest is susceptible to a pathotype and the differential is not susceptible  #
      # to that pathotype, remove that differential/gene from the list of candidate genes                 #
      if(isTRUE(df_trim[df_trim$Cultivar == cult, c] == "s" && df_trim[r,c] == "r") == TRUE){
        tocut <- append(tocut, df_trim[r,2])
      }}}
  
  ##### Remove varieties containing impossible genes, stored in tocut, from df_trim
  
  # Trim tocut to only unique values and remove those from df_trim
  tocut <- unique(tocut)
  df_trim <- df_trim[!(df_trim$Cultivar %in% tocut),]
  
  #### Sort df, df_raw, df_trim and df_trim_raw by matches and prepare for output
  
  #sort df_trim by matches
  df_trim <- df_trim[order(df_trim$Matches, decreasing = TRUE),]
  
}

################################################################################
############################### SERVER #########################################
################################################################################
# Define server logic to summarize and view selected dataset ----
server <- function(input, output, session) {
  
  # Reactive function to return the selected dataset  ----
  df <- reactive({
    req(input$df)
    read.csv(input$df$datapath)
  })
  
  df_raw <- reactive({
    req(input$df_raw)
    read.csv(input$df_raw$datapath)
  })
  
  # Update the cultivar list based on the uploaded data ----
  observe({
    req(df())
    updateSelectInput(session, "cult", choices = unique(df()$Cultivar))
  })
  
  
  # Call the appropriate function to return the selected dataset ----
  output$view <- renderTable({
    req(input$format)
    req(input$cult)
    req(input$obs)
    
    if (input$format == "Gene postulations (R/S)") {
      rsscores_trimmed(df_raw(), df(), input$cult)[1:input$obs,]
    } else if (input$format == "Gene postulations (Raw scores)") {
      rawscores_trimmed(df_raw(), df(), input$cult)[1:input$obs,]
    } else if (input$format == "Full sorted dataset (R/S)") {
      rsscores_sorted(df_raw(), df(), input$cult)[1:input$obs,]
    } else if (input$format == "Full sorted dataset (Raw scores)") {
      rawscores_sorted(df_raw(), df(), input$cult)[1:input$obs,]
    }
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0(input$cult, "_", input$format, ".csv")
    },
    content = function(file) {
      
      write.csv(
        switch(input$format,
               "Gene postulations (R/S)" = rsscores_trimmed(df_raw(), df(), input$cult)[1:input$obs,],
               "Gene postulations (Raw scores)" = rawscores_trimmed(df_raw(), df(), input$cult)[1:input$obs,],
               "Full sorted dataset (R/S)" = rsscores_sorted(df_raw(), df(), input$cult)[1:input$obs,],
               "Full sorted dataset (Raw scores)" = rawscores_sorted(df_raw(), df(), input$cult)[1:input$obs,]
        ),
        file,
        row.names = FALSE
      )
      
    }
  )
}