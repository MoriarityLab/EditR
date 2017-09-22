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
      Check the QC plot of signal / noise post filtering to pick values -- mouse over the plot for index numbers.
      Note: aggressive trimming may lead to inaccurate P-values, make sure to trim to the region that represents the majority of your sample."),
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
                         p("This plot shows the percent area of the signal for each base (ACGT) at each position along
                           the guide. Bases that are significantly different from the noise are colored in, color coded
                           relative to their percent area. Most positions of the guide only have one base colored in, as there
                           this is due to only one peak being present at that position. In the case of editing, there are more 
                           than one base colored in at a position."),
                         h3("Base info"),
                         p("You can use the information in the following table to better understand the values given for the percent
                           area in the 'Table plot'. The average percent signal shows you the average percent area (area of signal / all
                           other base areas), and (100 - average percent signal) shows you how much noise is present during the reading 
                           of that base. The amount of noise present in each base is also reflected in the 'Critical percent value', 
                           which is the critical value for the percent area for that base -- any percent area value above that would 
                           be called significantly different from the noise, and therefore editing has occurred. This critical value is 
                           calculated based on the P-value cutoff specified (default is 0.01). The model mu is the mu term for the zero
                            adjusted gamma distribution to model the noise -- a higher value of mu means that the predicted editing value
                            is less accurate. Filliben's correlation is how well the noise is modelled by a zero adjusted gamma 
                           distribution. Lower values of Filliben's correlation means less confidence in EditR for predicting editing,
                           your values should be above 0.90."),
                         tableOutput(outputId = "baseinfo.table"), 
                         h2("Quad plot"),
                         plotOutput(outputId = "editing.quad.plot"),
                         p("This plot shows the percent area of the background bases (the bases that 
                           were not the guide squence), the line denotes the crtical value cutoff based on the P-value
                           cutoff specified (default is 0.010.")
                ),                
                
                
                #### Report download tab
                tabPanel("Download Report",
                         h2("Download html report here"),
                         p("Click the button below to download a html report of the
                           predicted editing for the sequence that you uploaded. The report contains table output 
                          of the guide region, with the calculated probability of each base being from the noise at each position
                          and the base info table as well. Additionally, the editing table and the base information table are 
                          output in R syntax so you can copy them and bring them into your own R session.
Please be patient, as the report may take some time to download as it needs to be compiled."),
                         downloadButton('downloadReport')
                         )
                
                )
    

)
  
  )
)
