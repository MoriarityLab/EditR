shinyUI(
  
  fluidPage(

  titlePanel("EditR -- Version 2017-09-18_alpha"),
  
  sidebarPanel(
    p("You will need to supply your Sanger sequencing results file, as well as the guide sequence.
      Data are not stored on the server after the analysis is run."),
    checkboxInput("example", "Load Example Data", FALSE),
    
    conditionalPanel(condition = "!input.example",
                     p("Please upload your Sanger sequencing file here:"),
                     fileInput(inputId = "file", 
                               label = 'Upload .ab1 File'),
                     p("Please enter in your ~20bp guide sequence in 5'-3' orientation here. Should only contain: ACGT"),
                     textInput(inputId = 'guide',
                               label = 'Enter gRNA sequence'),
                     checkboxInput("guide.is.reverseComplement",
                                   "Guide sequence is reverse complement", FALSE)
                     ),
    
    p("Specify custom amounts to trim the input sanger sequence by.
      Check the QC plot of signal / noise post filtering to pick values -- mouse over the plot for index numbers"),
    numericInput(inputId ="trim5", label = "5' start", value = NA),
    numericInput(inputId ="trim3", label = "3' end", value = NA),
    p("Please enter your cutoff probability (p-value for significance), see instructions for details"),
      numericInput("pvalcutoff", "P-value cutoff:", 0.01, min = 0, max = 1)
  ),
  
  mainPanel(
    
    tabsetPanel(type = "tabs",
                
                #### Instructions tab
                #source("./instructions-content.R", local = TRUE)$value,
                # tabPanel("Instructions",
                # includeHTML("instructions.html")),
                
                tabPanel("Instructions",
                         includeMarkdown("instructions.md")),
                
                #### Data QC tab
                tabPanel("Data QC",
                         h3("Total peak area before filtering"),
                         plotOutput(outputId = "prefilter.totalarea"),
                         p("This plot shows the total peak area at each positon before QC filtering. 
                           The peak are for each base (ACGT) was summed to produce the total 
                           peak area. The shaded region denotes the guide RNA"),
                         
                         h3("Data QA: Signal and noise plot"),
                         plotlyOutput(outputId = "postfilter.signal.noise"),
                         textOutput(outputId = "trimmedrange"),
                         br(),
                         p(strong("Signal and noise plot:"), " This plot shows the peak area of the filtered dataset. This can tell
                           you if there are particular regions in your sequencing that have high
                           noise, and that it might affect the sensitivity of the predicted 
                           editing that we can make. The shaded region denotes the guide RNA region.
                           This plot shows the data postfiltering that will be used to generate the null
                           distribtution for our prediction method. The first 20 bases are removed as they
                           are generally poor quality. Then we set a filtering cutoff based on the mean of
                           the toal peak area, we remove positions where the total peak area is less than one 
                           tenth of the mean. For example, if the mean total peak area is 1000, we remove positions
                           that have less than 100 total peak area from our dataset. This mainly serves to trim the
                           low quality bases at the tail end of the sequencing run."),
                         
                         h3("Data QA: Percent noise peak area"),
                         plotOutput(outputId = "postfilter.noise.perc"),
                         p("This plot shows the amount of noise at a base position as the percent of the total peak area.
                            The shaded region is where the gRNA matches, which might have a peak due to base editing."),
                         h3("Guide RNA region Chromatogram"),
                         plotOutput(outputId = "chromatogram")
                ),
                
                #### Predicted editing tab
                tabPanel("Predicted Editing", 
                         h2("Table plot"),
                         plotOutput(outputId = "editing.table.plot"),
                         h2("Plot of where editing has likely ocurred"),
                         p("Line denotes significance cutoff."),
                         plotOutput(outputId = "editing.quad.plot")
                ),                
                
                
                #### Report download tab
                tabPanel("Download Report",
                         h2("Download html report here"),
                         p("Click the button below to download a html report of the
                           predicted editing for the sequence that you uploaded. The
                           report may take some time to generate as all the calculations
                           need to be rerun. Please be patient."),
                         downloadButton('downloadReport')
                         )
                
                )
    

)
  
  )
)
