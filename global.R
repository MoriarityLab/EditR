# version
editrversion = "1.0.10"

# loading all the packages into the global environment

library(shiny)
library(Biostrings)
library(sangerseqR)
library(gamlss)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(rmarkdown)
library(plotly)
library(yaml)


## Defining functions



GetGuideMatch <- function(guide, input.seq)
  # this function takes in the guide sequence, and the primary seq from a makeBaseCalls object
  # guide should be a DNAString object -- DNAString(guideinput)
  # input.seq = input.basecalls@primarySeq
  # returns guide.coord, a list with $match (forward / reverse), $start, and $end
  # and finds where the guide matches
  # it looks for both a forward and reverse match
  # in deciding if it's forward / reverse, it takes the match that is the longest
  # This function tries to extend the guide coordinates in case there is a partial match -- this might break
  # match to guide
  {

  
  # alignments
  align.f  <- pairwiseAlignment(pattern=guide, subject=input.seq, type="overlap")

  guide.coord <- list(match = "forward", start = align.f@subject@range@start,
                      end = align.f@subject@range@start + align.f@subject@range@width-1)
  
  # check to see if a match was actually found
  if(align.f@pattern@range@width <= 1){
    stop("Failed to find a forward match of the gRNA -- should it be a reverse complement?")
    return(NA)
  }
  

  # check to see if the matched region == the guide length
  is.match.guide.length <- (guide.coord$end - guide.coord$start + 1) == length(guide)
  
  # otherwise figure out what part of the guide matched to where, and update the guide region to include the unmatched areas. 
  if(is.match.guide.length != TRUE){
    
    # find the diff from beginning
    startdiff <- 1 - align.f@pattern@range@start
    # find the diff from the end
    enddiff <- 20 - (align.f@pattern@range@start + align.f@pattern@range@width-1)
    # add these numbers to the guide.cord
    guide.coord$start <- guide.coord$start + startdiff
    guide.coord$end <- guide.coord$end + enddiff
    
    # check to make sure that I didn't do something bad
    # if(guide.coord$start < 1 | guide.coord$end > length(guide)){
    #   stop("Something went wrong in matching the guide to the sequence, partial match of the guide found.
    #          Error in extending the match to the full guide length.")
    # }
    
    
  }
  
  return(guide.coord)
}


CreateSangs <- function(peakAmp, basecalls){
  # this function takes the peak amp matrix, and the input basecalls, and creates an 
  # unfiltered sanger sequencing dataframe (that I like to call sangs)
  
  
  sangs <- as.data.frame(peakAmp)
  names(sangs) <- c("A.area","C.area","G.area","T.area")
  
  # adding total area, and then calculating the percent area of each base
  sangs %<>% 
    mutate(Tot.area = A.area + C.area + G.area + T.area,
           A.perc = 100*A.area / Tot.area,
           C.perc = 100*C.area / Tot.area,
           G.perc = 100*G.area / Tot.area,
           T.perc = 100*T.area / Tot.area) 
  
  # adding on base calls
  sangs$base.call <- strsplit(x = toString(basecalls@primarySeq), split = "") %>%
    unlist
  
  # adding an index
  sangs$index <- seq_along(sangs$base.call)
  return(sangs)
}


GetNullDistModel <- function(sangs.filt, guide.coord)
  # this function takes the sangs.filt dataframe, gathers values for the null distribution
  # for each base, and then fits a zero-adjusted gamma distribution to these values
  # it returns a list of data.frames, each with the parameters of a null distribution 
  # for each base
  {

  
  sangs.filt %<>% filter(!(index %in% (guide.coord$start:guide.coord$end)) )
  
  nvals <- list()
  nvals$t <- sangs.filt %>% filter(base.call != "T") %>% select(T.perc) %>% unlist()
  nvals$c <- sangs.filt %>% filter(base.call != "C") %>% select(C.perc) %>% unlist()
  nvals$g <- sangs.filt %>% filter(base.call != "G") %>% select(G.perc) %>% unlist()
  nvals$a <- sangs.filt %>% filter(base.call != "A") %>% select(A.perc) %>% unlist()
  
  # Updated 3.26.19 to account for ultra clean sequencing
  replacement_zaga = c(rep(0, 989), 0.00998720389310502, 0.00998813447664401,0.009992887520785,
                       0.00999585366068316, 0.00999623914632598, 0.00999799013526835, 0.010001499423723,
                       0.0100030237039207, 0.0100045782875701, 0.0100048452355807, 0.0100049548867042)
  
  n.models <- lapply(nvals, FUN = function(x){
    set.seed(1)
    if((unique(x)[1] == 0 & length(unique(x)) == 1) |
       (unique(x)[1] == 0 & length(unique(x)) == 2 & table(x)[2] == 1))
    {x = replacement_zaga; message("Replacement vector used for low noise.")} # add noise if all 0s, or all 0s and one other value.
    tryCatch(gamlss((x)~1, family = ZAGA), error=function(e) # Progressively step up the mu.start if it fails
      tryCatch(gamlss((x)~1, family = ZAGA, mu.start = 1), error=function(e) 
        tryCatch(gamlss((x)~1, family = ZAGA, mu.start = 2), error=function(e) 
          tryCatch(gamlss((x)~1, family = ZAGA, mu.start = 3), error=function(e) # additional step added.
            gamlss((x)~1, family = ZAGA, mu.start = mean(x))
          )
        )
      )
    )
    # throws errors when a completely 0 vector
  })
  
  null.m.params <- lapply(n.models, FUN = function(x){
    mu <- exp(x$mu.coefficients[[1]])
    sigma <- exp(x$sigma.coefficients[[1]])
    nu.logit <- x$nu.coefficients[[1]]
    nu <- exp(nu.logit)/(1+exp(nu.logit))
    fillibens <-cor(as.data.frame(qqnorm(x$residuals, plot = FALSE)))[1,2]
    
    return(data.frame(mu= mu, sigma = sigma, nu = nu, fillibens = fillibens))
  })
  
  return(null.m.params)
}


CreateEditingDF <- function(guide.coord, guide, sangs, null.m.params){
  # this function creates the editing.df that contains the prob that the peak areas 
  # are not part of the noise distibution
  
  guide.df <- sangs[guide.coord$start:guide.coord$end,]
  
  # adding the guide sequence onto the data.frame
  guide.df$guide.seq <- guide %>% toString() %>% strsplit(. , "") %>% unlist
  
  
  # usage (params = null.m.params$t, perc = guide.df$T.perc)
  calcBaseProb <- function(params, perc){
    outprobs <- pZAGA(q = perc,
                      mu = params$mu,
                      sigma = params$sigma,
                      nu = params$nu,
                      lower.tail = FALSE)
  }
  
  guide.df$T.pval <- calcBaseProb(params = null.m.params$t, perc = guide.df$T.perc)
  guide.df$C.pval <- calcBaseProb(params = null.m.params$c, perc = guide.df$C.perc)
  guide.df$G.pval <- calcBaseProb(params = null.m.params$g, perc = guide.df$G.perc)
  guide.df$A.pval <- calcBaseProb(params = null.m.params$a, perc = guide.df$A.perc)
  
  guide.df$guide.position <- seq_along(guide.df$guide.seq)
  
  editing.df <- guide.df  
  return(editing.df)
}


makeEditingBarPlot <- function(edit.long, null.m.params, base, pval, editing.df){
  # makeEditingBarPlot(edit.filt = filter(editing.df, guide.seq != "T"),
  # null.m.params = null.m.params$t, base = "T, pval = p.val.cutoff)
  
  
  cutoff.line <- qZAGA(p = pval, mu = null.m.params$mu,
                       sigma = null.m.params$sigma, 
                       nu = null.m.params$nu, lower.tail = FALSE)
  
  edit.filt <- edit.long %>% filter(guide.seq != base, focal.base == base)
  
  edit.filt %>% filter(measure == "perc") %>% 
    ggplot(aes(x = guide.position, y = value)) +
    geom_bar(stat = "identity") +
    scale_y_continuous(breaks = seq(0, 100, 10)) +
    coord_cartesian(ylim = c(0,100)) + 
    scale_x_continuous(breaks = seq_along(editing.df$guide.seq)) + 
    annotate(geom = "segment", x = 0, xend = length(editing.df$guide.seq),
             y = cutoff.line, yend = cutoff.line) +
    labs(title = paste0("Percent ", base, " editing"),
         x = "Base postion in guide RNA",
         y = "Percent Total peak area") + 
    theme(axis.text.x = element_text(size = 8))
  
}

