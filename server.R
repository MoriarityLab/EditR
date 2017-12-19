## having example data

# the example file was orginally sequenced from Moriarity Lab work,
  # original file was 3a_PVT1For
examplefile <- "example.ab1"
exampleguide <- "CACTGGAATGACACACGCCC"

shinyServer(
  
  function(input, output) {
    
    
    # reading in file
    input.seqReactive <- reactive({
      # if else things to handle example data loading
      if(input$example) {
        return(readsangerseq(examplefile))
      } else if(!is.null(input$file)) {
        return(readsangerseq(input$file$datapath))
      } else return(validate(
        need(input$file, "Please upload your sanger sequence file")
      ))
      
      
    })
    
    # making basecalls
    
    input.basecallsReactive <- reactive({
      makeBaseCalls(input.seqReactive())
    })
    
    # getting the guide sequence
    guideReactive <- reactive({
      
      
      if(input$example) {
        guide <- (DNAString(exampleguide))
        }
      else{ 
        validate(
          need(input$guide != "", "Please enter a guide RNA sequence")
        )
        
        guide <- (DNAString(input$guide))
        
        if(input$guide.is.reverseComplement){
          guide <- reverseComplement(guide)
        }
      } 
      
      return(guide)
    })
    
    ## pvalue cutoff
    
    p.val.Reactive <- reactive({
      validate(
        need(input$pvalcutoff != "", "Please enter a p-value cutoff")
      )
      
      input$pvalcutoff
    })
    
    
    
    sangsReactive <- reactive({
      input.seq <- input.seqReactive()
      input.basecalls <- input.basecallsReactive()
      input.peakampmatrix <- peakAmpMatrix(input.basecalls)
      
      
      ### creating a sanger sequencing data.frame
      sangs <- CreateSangs(input.peakampmatrix, input.basecalls)
    
      return(sangs)
    })

    
    sangs.filtReactive <- reactive({
      sangs <- sangsReactive()
      
      if(is.na(input$trim3) & is.na(input$trim5)){
        # Remove 5' poor sequencing due to poor primer binding
        sangs.filt <- sangs %>% filter(index > 20)
        # removing crappy end
        peakTotAreaCutoff <- mean(sangs.filt$Tot.area)/10
        sangs.filt %<>% filter(Tot.area > peakTotAreaCutoff)
        return(sangs.filt)
      }else{
        validate(need(!is.na(input$trim3), label = "A value to trim the 3' end "))
        validate(need(!is.na(input$trim5), label = "A value to trim the 5' end "))

        sangs.filt <- sangs[input$trim5:input$trim3, ]
        return(sangs.filt)
      }
      

    })
    
    
    output$trimmedrange <- renderText({
      
      if(!is.na(input$trim5) & !is.na(input$trim3)){
        paste("Trimmed input, starting at ", input$trim5, " and ending at ", input$trim3)
      }

    })
    

    # finding the guide coordinates
    guide.coordReactive <- reactive({
      # the guide coordinates are relative to the index of the sequence that came out of the 
      # makeBaseCalls function
      # the guide is matched to the sequence of the filtered sanger sequencing, so that it doesn't
      # match to any of the low quality regions. 
      
      sangs.filt <- sangs.filtReactive()
      
      # getting the sequence to match to
      filt.sequence <- sangs.filt$base.call %>% paste(collapse = "") %>% DNAString()
      
      # this function finds where the guide matches 
      guide.match <- GetGuideMatch(guideReactive(), filt.sequence)
      
      # Finding the index values
      guide.coord <- list(start = sangs.filt[guide.match$start, "index"],
                          end = sangs.filt[guide.match$end, "index"])
      
      
      return(guide.coord)
      
    })
    
    


  
    #######################
    ## Model fitting work
    
    
    nullparams.Reactive <- reactive({
      sangs.filt <- sangs.filtReactive()
      guide.coord <- guide.coordReactive()
      # getting params for the different null models
      null.m.params <- GetNullDistModel(sangs.filt, guide.coord)
    })
    
    
    # Editing reactive
    # this produces a dataframe of the guide region with the probabilities that the nonprimary base shows evidence of editing
    
    editing.Reactive <- reactive({
      # grabbing the reactive stuff
      sangs.filt <- sangs.filtReactive()
      sangs <- sangsReactive()
      guide.coord <- guide.coordReactive()
      guide <- guideReactive()
      null.m.params <- nullparams.Reactive()
      
      
      editing.df <- CreateEditingDF(guide.coord, guide, sangs, null.m.params)
      
      return(editing.df)
    })
    
    
    ### Making QC plots
    
    output$prefilter.totalarea <- renderPlot({
      guide.coord <- guide.coordReactive()
      sangs <- sangsReactive()
      ggplot(sangs, aes(x = index, y = Tot.area)) +
        geom_rect(xmin = guide.coord$start, xmax = guide.coord$end, ymin = 0, ymax = Inf, fill = "lightgrey") +
        geom_line() +
        labs(x = "Base position",
             y = "Total peak area at base position (A + T + C + G)",
             title = "Unfiltered data, total peak area")
    })
    
    output$postfilter.signal.noise <- renderPlotly({
      guide.coord <- guide.coordReactive()
      sangs.filt <- sangs.filtReactive()
      sangs <- sangsReactive()
      
      ## doing some data rearrangement for plotting
      
      catch <- sangs.filt %>% dplyr::select(A.area:T.area, A.perc:T.perc, base.call, index)
      catch %<>% gather( key = base, value = value,
                         A.area:T.area, A.perc:T.perc) %>%
        separate(col = base, into = c("base", "measure"))
      
      noise <- catch %>% 
        group_by(index, measure) %>% 
        filter(base != base.call) %>%
        summarize(noise = sum(value))
      signal <- catch %>%
        group_by(index, measure) %>%
        filter(base == base.call) %>% 
        summarize(signal = sum(value))
      
      summary.df <- left_join(noise,signal) %>%
        gather(key = type, value = value, noise, signal)
      
      sangs.plot <- sangs.filt %>% dplyr::select(index, base.call) 
      sangs.plot %<>% left_join(summary.df)
      
      ## now acutally making the plot
      
      sangs.plot$type <- ordered(sangs.plot$type, levels = c("signal", "noise"))
      
      max.height <- filter(sangs.plot, measure == "area") %>% select(value) %>% max
      
      guide.region <- data.frame(index = c(guide.coord$start, guide.coord$end,
                                           guide.coord$end, guide.coord$start), 
                                 value = c(0, 0, max.height, max.height))
      
      p <- sangs.plot %>% filter(measure == "area") %>%
        ggplot(aes(x = index, y = value)) +
        geom_polygon(data = guide.region, fill = "lightgrey") +
        geom_line(aes(color = type)) +
        scale_color_manual(values = c("#5e3c99", "#e66101")) +  
        labs(title = "Filtered data signal and noise total area",
             x = "Position",
             y = "Peak Area") + 
        theme(legend.title = element_blank())
      ggplotly(p)
      
    })
    
    output$postfilter.noise.perc <- renderPlot({
      guide.coord <- guide.coordReactive()
      sangs.filt <- sangs.filtReactive()
      
      guide.coord <- guide.coordReactive()
      sangs.filt <- sangs.filtReactive()
      
      ## doing some data rearrangement for plotting
      
      catch <- sangs.filt %>% dplyr::select(A.area:T.area, A.perc:T.perc, base.call, index)
      catch %<>% gather( key = base, value = value,
                         A.area:T.area, A.perc:T.perc) %>%
        separate(col = base, into = c("base", "measure"))
      
      # splitting the catch dataframe into either signal and noise, and calculating the 
      # total noise area or total noise percent
      noise <- catch %>% 
        group_by(index, measure) %>% 
        filter(base != base.call) %>%
        summarize(noise = sum(value))
      signal <- catch %>%
        group_by(index, measure) %>%
        filter(base == base.call) %>% 
        summarize(signal = sum(value))
      
      signal.noise.df <- left_join(noise,signal) %>%
        gather(key = type, value = value, noise, signal)
      
      # making the plotting df
      sangs.plot <- sangs.filt %>% dplyr::select(index, base.call) 
      sangs.plot %<>% left_join(signal.noise.df)
      sangs.plot$type <- ordered(sangs.plot$type, levels = c("signal", "noise"))
      
      sangs.plot %>% filter(measure == "perc", type == "noise") %>%
        ggplot(aes(x = index, y = value)) +
        geom_area(aes(fill = type)) +
        scale_fill_manual(values = c("#e66101")) + 
        annotate("rect", xmin=guide.coord$start, xmax=guide.coord$end,
                 ymin=0, ymax=Inf, alpha=1/5, fill="black") +
        guides(fill = FALSE) + 
        labs(title = "Percent peak area noise",
             x = "Position",
             y = "Percent peak area noise")
    })
    
    
    output$chromatogram <- renderPlot({
      guide.coord <- guide.coordReactive()
      input.basecalls <- input.basecallsReactive()
      
      chromatogram(obj = input.basecalls, 
                   trim5 = (guide.coord$start-1), trim3 = length(input.basecalls@primarySeq) - guide.coord$end,
                   width = (guide.coord$end - guide.coord$start + 1))
    })
    
    
    base.infoReactive <- reactive({
      # this reactive is collecting information on each base:
      # - average percent signal
      # critical value for each base for the cutoff of significance
      
      sangs.filt <- sangs.filtReactive()
      null.m.params <- nullparams.Reactive()
      p.val.cutoff <- p.val.Reactive()
      
      # finding the average percent signal for each base
      avg.base <- sangs.filt %>% gather(key = focal.base, value = value, 
                                A.area:T.area, A.perc:T.perc) %>%
        separate(col = focal.base, into = c("focal.base", "measure")) %>% 
        spread(key = measure, value = value) %>% 
        filter(base.call == focal.base) %>% 
        group_by(focal.base) %>% 
        summarize(avg.percsignal = mean(perc),
                  avg.areasignal = mean(area))
      
      # getting the critical value for each base
      crit.vals <- c(
        a = qZAGA(p = p.val.cutoff, mu = null.m.params$a$mu,
            sigma = null.m.params$a$sigma,
            nu = null.m.params$a$nu,
            lower.tail = FALSE),
        c = qZAGA(p = p.val.cutoff, mu = null.m.params$c$mu,
                  sigma = null.m.params$c$sigma,
                  nu = null.m.params$c$nu,
                  lower.tail = FALSE),
        g = qZAGA(p = p.val.cutoff, mu = null.m.params$g$mu,
                  sigma = null.m.params$g$sigma,
                  nu = null.m.params$g$nu,
                  lower.tail = FALSE),
        t = qZAGA(p = p.val.cutoff, mu = null.m.params$t$mu,
                  sigma = null.m.params$t$sigma,
                  nu = null.m.params$t$nu,
                  lower.tail = FALSE)
        )
      # getting fillibens
      fil <- lapply(null.m.params, FUN = function(x){x$fillibens})
      filvec <- c(a = fil$a, c = fil$c, g = fil$g, t = fil$t)
      
      # getting mu -- a measure of dispersion
      mul <- lapply(null.m.params, FUN = function(x){x$mu})
      mulvec <- c(a = mul$a, c = mul$c, g = mul$g, t = mul$t)
      
      return(data.frame(avg.base, crit.perc.area = crit.vals, mu = mulvec, fillibens = filvec))
    })
    
    
    output$baseinfo.table <- renderTable({
      base.info <- base.infoReactive()
      
      temp <- base.info
      names(temp) <- c("Base", "Average percent signal",
                       "Average peak area",  "Critical percent value",
                       "model mu",  "Fillibens correlation")
      row.names(temp) <- NULL
      
      return(temp)
    })
    
    ##############################################
    output$editing.table.plot <- renderPlot({
      editing.df <- editing.Reactive()
      null.m.params <- nullparams.Reactive()
      p.val.cutoff <- p.val.Reactive()
      
      
      edit.long <- editing.df %>% gather(key = focal.base, value = value, 
                                         A.area:T.area, A.perc:T.perc, T.pval:A.pval) %>%
        separate(col = focal.base, into = c("focal.base", "measure"))
      
      edit.spread <- edit.long %>% spread(key = measure, value = value)
      
      edit.color <- edit.spread %>% filter(pval < p.val.cutoff)
      
      edit.spread %>% 
        ggplot(aes(x = as.factor(index), y = focal.base)) + 
        geom_tile(data = edit.color, aes(fill = perc)) + 
        geom_text(aes(label = round(perc, 0)), angle = 0, size = 4) +   
        guides(fill = FALSE) + 
        scale_fill_continuous(low = "#f7a8a8", high = "#9acdee") + 
        scale_x_discrete(position = "top", labels = editing.df$guide.seq) + 
        labs(x = "Guide sequence", y = NULL) + 
        theme(axis.ticks = element_blank())
      
    })
    
###########################################
    output$editing.quad.plot <- renderPlot({
      editing.df <- editing.Reactive()
      null.m.params <- nullparams.Reactive()
      p.val.cutoff <- p.val.Reactive()
      
      edit.long <- editing.df %>% gather(key = focal.base, value = value, 
                                         A.area:T.area, A.perc:T.perc, T.pval:A.pval) %>%
        separate(col = focal.base, into = c("focal.base", "measure"))
      
      p.a <- makeEditingBarPlot(edit.long = edit.long, null.m.params = null.m.params$a,
                                base = "A", pval = p.val.cutoff, editing.df)
      p.c <- makeEditingBarPlot(edit.long, null.m.params$c,
                                base = "C", pval = p.val.cutoff, editing.df)
      p.g <- makeEditingBarPlot(edit.long, null.m.params$g,
                                base = "G", pval = p.val.cutoff, editing.df)
      p.t <- makeEditingBarPlot(edit.long, null.m.params$t,
                                base = "T", pval = p.val.cutoff, editing.df)
      
      grid.arrange(p.a, p.c, p.g, p.t)
    })



##########################33
    output$editingtable <- renderTable(digits = -3, expr = { # maybe need to use renderMarkdown?
      editing.df <- editing.Reactive()
      p.val.cutoff <- p.val.Reactive()
      
      edit.long <- editing.df %>% gather(key = focal.base, value = value, 
                                         A.area:T.area, A.perc:T.perc, T.pval:A.pval) %>%
        separate(col = focal.base, into = c("focal.base", "measure"))
      
      edit.spread <- edit.long %>% spread(key = measure, value = value)
    
      
      editingtable <- edit.spread %>% arrange(guide.position) %>% select(index,guide.position, guide.seq, base.call, focal.base, perc, pval) %>% 
        mutate(index = as.character(index), 
               guide.position = as.character(guide.position),
               perc = format(perc, digits = 2),
               signficant = ifelse(pval < p.val.cutoff, yes = "*", no = ""))
      
      names(editingtable) <-  c("Sanger position", "Guide position","Guide sequence", "Sanger base call",
                         "Focal base", "Focal base peak area", "p value", "")
      
      return(editingtable)
    })
    
    output$modelfits <- renderTable({
      null.m.params <- nullparams.Reactive()
      
      temp <- data.frame(`Base` = c("A", "C", "G", "T"),
                 `Fillibens Correlation` = c(null.m.params$a$fillibens, 
                                          null.m.params$g$fillibens, 
                                          null.m.params$c$fillibens, 
                                          null.m.params$t$fillibens))
      
    names(temp) <- c("Base", "Fillibens Correlation")
    return(temp)
    })
    
    
#################################3
    
    output$downloadReport <- downloadHandler(
      # For PDF output, change this to "report.pdf"
      filename = "report.html",
      content = function(file) {
        # Copy the report file to a temporary directory before processing it, in
        # case we don't have write permissions to the current working dir (which
        # can happen when deployed).
        tempReport <- file.path(tempdir(), "report.Rmd")
        file.copy("report.Rmd", tempReport, overwrite = TRUE)
        
        # Set up parameters to pass to Rmd document
         # pulling in the reactive things to generate the report with:
        
        params <- list(null.m.params = nullparams.Reactive(),
                       guide = guideReactive(),
                       guide.coord = guide.coordReactive(),
                       editing.df = editing.Reactive(),
                       p.val.cutoff = p.val.Reactive(),
                       input.seq = input.seqReactive(),
                       sangs = sangsReactive(),
                       sangs.filt = sangs.filtReactive(),
                       base.info = base.infoReactive(),
                       editrversion = editrversion)
        
        # Knit the document, passing in the `params` list, and eval it in a
        # child of the global environment (this isolates the code in the document
        # from the code in this app).
        rmarkdown::render(tempReport, output_file = file,
                          params = params,
                          envir = new.env(parent = globalenv())
        )
      }
    )

})
