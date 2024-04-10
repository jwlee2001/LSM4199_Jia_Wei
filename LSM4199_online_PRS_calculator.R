## breast cancer PRS tool ##

library(shinyjs)
library(shiny)
library(plotly)
library(survival)

## load dataset
load("brca_train_27_genes_scaled.RData")

## define UI

ui <- fluidPage(
  titlePanel("Breast Cancer Polygenic Risk Score Calculator"),
  h6(div(HTML("Use the sliding bars below to indicate the scaled TPM values for the 27 genes"))),
  useShinyjs(),
  #create client side input form
  sidebarLayout(
    sidebarPanel(
                   sliderInput(inputId = 'ENSG00000126353.3', label = 'C-C Motif Chemokine Receptor 7 (CCR7)', min=-5, max=20,value=0, step=0.1, round=0),
                   sliderInput(inputId = 'ENSG00000171951.5', label = 'Secretogranin II (SCG2)', min=-5, max=20,value=0, step=0.1, round=0),
                   sliderInput(inputId = 'ENSG00000124575.7', label = 'H1.3 Linker Histone, Cluster Member (H1-3)', min=-5, max=20,value=0, step=0.1, round=0),
                   sliderInput(inputId = 'ENSG00000185753.13', label = 'Chromosome X Open Reading Frame 38 (CXorf38)', min=-5, max=20,value=0, step=0.1, round=0),
                   sliderInput(inputId = 'ENSG00000160678.12', label = 'S100 Calcium Binding Protein A1 (S100A1)', min=-5, max=20,value=0, step=0.1, round=0),
                   sliderInput(inputId = 'ENSG00000128394.17', label = 'Apolipoprotein B MRNA Editing Enzyme Catalytic Subunit 3F (APOBEC3F)', min=-5, max=20,value=0, step=0.1, round=0),
                   sliderInput(inputId = 'ENSG00000185527.12', label = 'Phosphodiesterase 6G (PDE6G)', min=-5, max=20,value=0, step=0.1, round=0),
                   sliderInput(inputId = 'ENSG00000091106.19', label = 'NLR Family CARD Domain Containing 4 (NLRC4)', min=-5, max=20,value=0, step=0.1, round=0),
                   sliderInput(inputId = 'ENSG00000112339.15', label = 'HBS1 Like Translational GTPase (HBS1L)', min=-5, max=20,value=0, step=0.1, round=0),
                   sliderInput(inputId = 'ENSG00000095739.11', label = 'BMP And Activin Membrane Bound Inhibitor (BAMBI)', min=-5, max=20,value=0, step=0.1, round=0),
                   sliderInput(inputId = 'ENSG00000140374.16', label = 'Electron Transfer Flavoprotein Subunit Alpha (ETFA)', min=-5, max=20,value=0, step=0.1, round=0),
                   sliderInput(inputId = 'ENSG00000261068.2', label = 'Clone-based (Vega) Gene', min=-5, max=20,value=0, step=0.1, round=0),
                   sliderInput(inputId = 'ENSG00000148926.10', label = 'Adrenomedullin (ADM)', min=-5, max=20,value=0, step=0.1, round=0),
                   sliderInput(inputId = 'ENSG00000272668.2', label = 'Novel Transcript (LOC107985216)', min=-5, max=20,value=0, step=0.1, round=0),
                   sliderInput(inputId = 'ENSG00000102078.16', label = 'Solute Carrier Family 25 Member 14 (SLC25A14)', min=-5, max=20,value=0, step=0.1, round=0),
                   sliderInput(inputId = 'ENSG00000243679.1', label = 'Novel Pseudogene (ENSG00000243679)', min=-5, max=20,value=0, step=0.1, round=0),
                   sliderInput(inputId = 'ENSG00000104267.10', label = 'Carbonic Anhydrase 2 (CA2)', min=-5, max=20,value=0, step=0.1, round=0),
                   sliderInput(inputId = 'ENSG00000183807.8', label = 'Family With Sequence Similarity 162 Member B (FAM162B)', min=-5, max=20,value=0, step=0.1, round=0),
                   sliderInput(inputId = 'ENSG00000226038.5', label = 'Peptidylprolyl Isomerase A Pseudogene 21', min=-5, max=20,value=0, step=0.1, round=0),
                   sliderInput(inputId = 'ENSG00000254635.6', label = 'WAC Antisense RNA 1 (Head To Head)', min=-5, max=20,value=0, step=0.1, round=0),
                   sliderInput(inputId = 'ENSG00000174327.7', label = 'Solute Carrier Family 16 Member 13 (SLC16A13)', min=-5, max=20,value=0, step=0.1, round=0),
                   sliderInput(inputId = 'ENSG00000089723.10', label = 'OTU Deubiquitinase, Ubiquitin Aldehyde Binding 2 (OTUB2)', min=-5, max=20,value=0, step=0.1, round=0),
                   sliderInput(inputId = 'ENSG00000152465.18', label = 'N-Myristoyltransferase 2 (NMT2)', min=-5, max=20,value=0, step=0.1, round=0),
                   sliderInput(inputId = 'ENSG00000163666.10', label = 'HESX Homeobox 1 (HESX1)', min=-5, max=20,value=0, step=0.1, round=0),
                   sliderInput(inputId = 'ENSG00000253549.6', label = 'CA3 Antisense RNA 1', min=-5, max=20,value=0, step=0.1, round=0),
                   sliderInput(inputId = 'ENSG00000136999.5', label = 'Cellular Communication Network Factor 3 (CCN3)', min=-5, max=20,value=0, step=0.1, round=0),
                   sliderInput(inputId = 'ENSG00000068024.17', label = 'Histone Deacetylase 4 (HDAC4)', min=-5, max=20,value=0, step=0.1, round=0),
                 ),
                 mainPanel(
                   # Output: Interactive histogram plot
                   plotlyOutput(outputId = "histogram"),
                   hr(),
                   # Display PRS value
                   textOutput("text"), 
                   # Display probability
                   textOutput("text2")
                 )
    )
  )

## define server logic

server <- function(input, output) {
  ## calculate PRS
  prs_calculate <- reactive({
    as.numeric(input$ENSG00000126353.3)*6.128737525 +
      as.numeric(input$ENSG00000171951.5)*1.370759151 +
      as.numeric(input$ENSG00000124575.7)*0.916218582 +
      as.numeric(input$ENSG00000185753.13)*-1.404059217 +
      as.numeric(input$ENSG00000160678.12)*1.503405398 +
      as.numeric(input$ENSG00000128394.17)*-1.520620105 +
      as.numeric(input$ENSG00000185527.12)*-1.680976092 +
      as.numeric(input$ENSG00000091106.19)*1.446239326 +
      as.numeric(input$ENSG00000112339.15)*0.851167948 +
      as.numeric(input$ENSG00000095739.11)*1.345985008 +
      as.numeric(input$ENSG00000140374.16)*0.962081293 +
      as.numeric(input$ENSG00000261068.2)*1.100483352 +
      as.numeric(input$ENSG00000148926.10)*1.346954533 +
      as.numeric(input$ENSG00000272668.2)*0.877209395 +
      as.numeric(input$ENSG00000102078.16)*1.688659281 +
      as.numeric(input$ENSG00000243679.1)*-1.36410439 +
      as.numeric(input$ENSG00000104267.10)*0.544201573 +
      as.numeric(input$ENSG00000183807.8)*0.830753561 +
      as.numeric(input$ENSG00000226038.5)*1.606417565 +
      as.numeric(input$ENSG00000254635.6)*1.217496297 +
      as.numeric(input$ENSG00000174327.7)*-1.183108945 +
      as.numeric(input$ENSG00000089723.10)*-1.282633421 +
      as.numeric(input$ENSG00000152465.18)*-0.964442315 +
      as.numeric(input$ENSG00000163666.10)*0.756811161 +
      as.numeric(input$ENSG00000253549.6)*-0.708942189 +
      as.numeric(input$ENSG00000136999.5)*1.066506179 +
      as.numeric(input$ENSG00000068024.17)*-1.021264492
  })
  
  # Render PRS
  output$text <- renderText({
    paste("The patient's PRS is:", round(prs_calculate(),digits=2))
  })
  
  # Function to calculate survival probability
  prob_calculate <- reactive({
    # Fit Cox proportional hazards model
    cox_model <- coxph(Surv(survival_time_years, censored_status_opp) ~ polygenic, data = brca_train_27_genes_scaled)
    # calculate PRS
    prs_temp = as.numeric(input$ENSG00000126353.3)*6.128737525 +
      as.numeric(input$ENSG00000171951.5)*1.370759151 +
      as.numeric(input$ENSG00000124575.7)*0.916218582 +
      as.numeric(input$ENSG00000185753.13)*-1.404059217 +
      as.numeric(input$ENSG00000160678.12)*1.503405398 +
      as.numeric(input$ENSG00000128394.17)*-1.520620105 +
      as.numeric(input$ENSG00000185527.12)*-1.680976092 +
      as.numeric(input$ENSG00000091106.19)*1.446239326 +
      as.numeric(input$ENSG00000112339.15)*0.851167948 +
      as.numeric(input$ENSG00000095739.11)*1.345985008 +
      as.numeric(input$ENSG00000140374.16)*0.962081293 +
      as.numeric(input$ENSG00000261068.2)*1.100483352 +
      as.numeric(input$ENSG00000148926.10)*1.346954533 +
      as.numeric(input$ENSG00000272668.2)*0.877209395 +
      as.numeric(input$ENSG00000102078.16)*1.688659281 +
      as.numeric(input$ENSG00000243679.1)*-1.36410439 +
      as.numeric(input$ENSG00000104267.10)*0.544201573 +
      as.numeric(input$ENSG00000183807.8)*0.830753561 +
      as.numeric(input$ENSG00000226038.5)*1.606417565 +
      as.numeric(input$ENSG00000254635.6)*1.217496297 +
      as.numeric(input$ENSG00000174327.7)*-1.183108945 +
      as.numeric(input$ENSG00000089723.10)*-1.282633421 +
      as.numeric(input$ENSG00000152465.18)*-0.964442315 +
      as.numeric(input$ENSG00000163666.10)*0.756811161 +
      as.numeric(input$ENSG00000253549.6)*-0.708942189 +
      as.numeric(input$ENSG00000136999.5)*1.066506179 +
      as.numeric(input$ENSG00000068024.17)*-1.021264492
    
    # new data
    new_data <- data.frame(
      survival_time_years = c(5),
      censored_status_opp = c(0),
      polygenic = c(prs_temp)
    )
    # Predict survival probabilities at the specified time
    predicted_survival <- predict(cox_model, newdata = new_data, type = "survival", se.fit = TRUE)
    return(predicted_survival$fit)
  })
  
  ## Render surv probability
  
  output$text2 <- renderText({
    paste("The patient's probability of 5-year survival is:", round(prob_calculate(),digits=2))
  })
   
  # Render histogram plot
  output$histogram <- renderPlotly({
    
    # Create plotly histogram with a vertical line at the summed value
    p <- plot_ly(x = brca_train_27_genes_scaled$polygenic, type = "histogram",name="Histogram of PRS") %>%
      layout(title = paste("Distribution of Polygenic Risk Scores"),
             xaxis = list(title = "Polygenic Risk Score"))
    
    # Add vertical line to indicate PRS
    p <- add_trace(p, x = c(prs_calculate(),prs_calculate()), y=c(0,150),type = "scatter", mode = "lines", name = "Patient's PRS", line = list(color = "red"))
    
    p
  })
}


## run the application
shinyApp(ui = ui, server = server)