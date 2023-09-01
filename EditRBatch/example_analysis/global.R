### global.R for EditR Batch
#############################################################################
# Copyright (C) 2023 Mitchell Kluesner (mkluesne@fredhutch.org)
#  
# This file is part of multiEditR (Multiple Edit Deconvolution by Inference of Traces in R)
# 
# Please only copy and/or distribute this script with proper citation of 
# multiEditR publication
#############################################################################

### Libraries
library(sangerseqR)
library(tidyverse)
library(magrittr)
library(plyr)
library(sangerseqR)
library(gamlss)
library(readr)
library(optparse)
library(readxl)
library(gridExtra)
library(DT)

### Data
bases = c("A", "C", "G", "T")
ACGT = bases

pre_trinucleotide_table = read_tsv("/Users/kluesner/Desktop/Research/EditR/multiEditR/pre_trinucleotide_table.tsv")

# # Load the GBM model
# load("/Users/kluesner/Desktop/Research/EditR/multiEditR/program/working_branch/gbm_model.rda")

phred_scores = data.frame(stringsAsFactors=FALSE,
                          phred = c("!", "“", "$", "%", "&", "‘", "(", ")", "*", "+", ",", "–",
                                    ".", "/", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9",
                                    ":", ";", "<", "=", ">", "?", "@", "A", "B", "C", "D", "E", "F",
                                    "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S",
                                    "T", "U", "V", "W", "X", "Y", "Z", "[", "\\", "]", "^", "_",
                                    "`", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l",
                                    "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y",
                                    "z", "{", "|", "}", "~"),
                          prob = c(1, 0.794328235, 0.501187234, 0.398107171, 0.316227766,
                                   0.251188643, 0.199526232, 0.158489319, 0.125892541, 0.1,
                                   0.079432824, 0.063095734, 0.050118723, 0.039810717, 0.031622777,
                                   0.025118864, 0.019952623, 0.015848932, 0.012589254, 0.01,
                                   0.007943282, 0.006309573, 0.005011872, 0.003981072, 0.003162278,
                                   0.002511886, 0.001995262, 0.001584893, 0.001258925, 0.001, 0.000794328,
                                   0.000630957, 0.000501187, 0.000398107, 0.000316228, 0.000251189,
                                   0.000199526, 0.000158489, 0.000125893, 1e-04, 7.94328e-05,
                                   6.30957e-05, 5.01187e-05, 3.98107e-05, 3.16228e-05, 2.51189e-05,
                                   1.99526e-05, 1.58489e-05, 1.25893e-05, 1e-05, 7.9433e-06,
                                   6.3096e-06, 5.0119e-06, 3.9811e-06, 3.1623e-06, 2.5119e-06, 1.9953e-06,
                                   1.5849e-06, 1.2589e-06, 1e-06, 7.943e-07, 6.31e-07, 5.012e-07,
                                   3.981e-07, 3.162e-07, 2.512e-07, 1.995e-07, 1.585e-07, 1.259e-07,
                                   1e-07, 7.94e-08, 6.31e-08, 5.01e-08, 3.98e-08, 3.16e-08,
                                   2.51e-08, 2e-08, 1.58e-08, 1.26e-08, 1e-08, 7.9e-09, 6.3e-09, 5e-09,
                                   4e-09, 3.2e-09, 2.5e-09, 2e-09, 1.6e-09, 1.3e-09, 1e-09, 8e-10,
                                   6e-10, 5e-10)
)

### Functions
# convert nucleotide to a factor
nucleotide_factor = function(x){factor(x, levels = c("A", "C", "G", "T"))}

# Function to analyze an individual .ab1 file
# Need to make sure the file ends in .ab1 for this function
make_ctrl_sanger_df = function(sanger_file){
  base_calls = makeBaseCalls(sanger_file)
  sanger_df = base_calls %>% peakAmpMatrix %>% data.frame()
  colnames(sanger_df) = c("A_area","C_area","G_area","T_area")
  sanger_df %<>% 
    mutate(., max_base = {apply(., 1, which.max) %>% bases[.]}) %<>%
    mutate(Tot.Area = A_area + C_area + G_area + T_area,
           A_perc = A_area / Tot.Area,
           C_perc = C_area / Tot.Area,
           G_perc = G_area / Tot.Area,
           T_perc = T_area / Tot.Area) %<>%
    mutate(base_call = strsplit(x = toString(base_calls@primarySeq), split = "") %>% unlist) %<>%
    mutate(index = 1:NROW(.)) %<>%
    mutate(max_base_height = {ifelse(max_base == "A", A_area,
                                     ifelse(max_base == "C", C_area,
                                            ifelse(max_base == "G", G_area,
                                                   ifelse(max_base == "T", T_area,NA))))})
}


### Makes a df for the edited sample that accounts for the misalignment of peaks
### Employs the secondary basecalls
### Not aligned for peaks that are shifted by a complete basecall
make_samp_sanger_df = function(samp_sanger, ctrl_seq){
  
  ### Align the phase of the primary and secondary basecalls in the sample sequence to that of the control
  ### This changes where and what the bases are called as
  ### Returns an object of class sangerseq
  phased_samp_sanger = samp_sanger %>%
    makeBaseCalls(.) %>%
    setAllelePhase(., ctrl_seq) # Does not require an actual chromatogram, so it is ammenable to just using a fasta file
  
  ### Return a data frame from the phased sample sangerseq object with the position of each 
  ### This method appears to return higher intensities for the noise, thus we're going to trust it more for detecting noise for modelling
  samp_peakAmpDF = phased_samp_sanger@peakPosMatrix %>%
    as_tibble() %>%
    dplyr::rename(A = V1, C = V2, G = V3, `T` = V4) %>%
    # na_if(., 0) %>%
    dplyr::mutate(primary_base_call = primarySeq(phased_samp_sanger) %>% as.character() %>% str_split(., pattern = "") %>% unlist(),
                  secondary_base_call = secondarySeq(phased_samp_sanger) %>% as.character() %>% str_split(., pattern = "") %>% unlist()) %>%
    dplyr::mutate(identical = {primary_base_call == secondary_base_call}) %>%
    dplyr::mutate(row_sum = as.integer(!is.na(A)) +
                    as.integer(!is.na(C)) +
                    as.integer(!is.na(G)) +
                    as.integer(!is.na(`T`))
    ) %>% 
    dplyr::mutate(peak_pos = {ifelse(grepl("A|C|G|T", primary_base_call),
                                     ifelse(primary_base_call == "A", A, 
                                            ifelse(primary_base_call == "C", C, 
                                                   ifelse(primary_base_call == "G", G, 
                                                          ifelse(primary_base_call == "T", `T`, "Error, stop!")))),
                                     pmin(A, C, G, `T`, na.rm = TRUE))}) %>%
    dplyr::mutate(A = ifelse(is.na(A), peak_pos, A) %>% as.numeric(),
                  C = ifelse(is.na(C), peak_pos, C) %>% as.numeric(),
                  G = ifelse(is.na(G), peak_pos, G) %>% as.numeric(),
                  `T` = ifelse(is.na(`T`), peak_pos, `T`) %>% as.numeric())
  
  ### Reformat df to be identical to output of make_ctrl_sanger_df()
  samp_peakAmpDF %<>%
    dplyr::mutate(A_area = samp_sanger@traceMatrix[samp_peakAmpDF$A, 1],
                  C_area = samp_sanger@traceMatrix[samp_peakAmpDF$C, 2],
                  G_area = samp_sanger@traceMatrix[samp_peakAmpDF$G, 3],
                  T_area = samp_sanger@traceMatrix[samp_peakAmpDF$`T`, 4]) %<>%
    dplyr::select(A_area:T_area, primary_base_call, secondary_base_call) %<>%
    mutate(., max_base = {apply(., 1, which.max) %>% bases[.]}) %<>%
    mutate(max_base = factor(max_base, levels = ACGT)) %<>%
    mutate(Tot.Area = A_area + C_area + G_area + T_area,
           A_perc = A_area / Tot.Area,
           C_perc = C_area / Tot.Area,
           G_perc = G_area / Tot.Area,
           T_perc = T_area / Tot.Area) %>%
    mutate(index = 1:NROW(Tot.Area)) %<>%
    mutate(max_base_height = {ifelse(max_base == "A", A_area,
                                     ifelse(max_base == "C", C_area,
                                            ifelse(max_base == "G", G_area,
                                                   ifelse(max_base == "T", T_area,NA))))})
  
  return(samp_peakAmpDF)
}

# The alignment index and subsequent filtering is currently off -- need to adapt to take the longest region of alignment, and then only filter our those sequences.
align_sanger_dfs = function(control_df, sample_df){
  
  ctrl_seq = control_df$max_base %>% paste0(., collapse = "") %>% DNAString()
  sample_seq = sample_df$max_base %>% paste0(., collapse = "") %>% DNAString()
  
  if(ctrl_seq@length > sample_seq@length){
    alignment = pairwiseAlignment(pattern = sample_seq, subject = ctrl_seq)
    control_df = control_df %>%
      mutate(align_index = 1:nrow(.)) %>%
      filter(align_index >= alignment@subject@range@start) %>%
      filter(align_index < (alignment@subject@range@start+alignment@subject@range@width))
    sample_df = sample_df %>%
      mutate(align_index = 1:nrow(.)) %>%
      filter(align_index >= alignment@pattern@range@start) %>%
      filter(align_index < alignment@pattern@range@start+alignment@pattern@range@width)
  } else
  {
    alignment = pairwiseAlignment(pattern = ctrl_seq, subject = sample_seq)
    control_df = control_df %>%
      mutate(align_index = 1:nrow(.)) %>%
      filter(align_index >= alignment@pattern@range@start) %>%
      filter(align_index < alignment@pattern@range@start+alignment@pattern@range@width)
    sample_df = sample_df %>%
      mutate(align_index = 1:nrow(.)) %>%
      filter(align_index >= alignment@subject@range@start) %>%
      filter(align_index < (alignment@subject@range@start+alignment@subject@range@width))
  }
  return(list(control_df, sample_df))
}

### Convert the abif to fastq
abif_to_fastq = function (seqname = "sample", path, trim = TRUE, cutoff = 1, 
                          min_seq_len = 20, offset = 33, recall = FALSE) 
{
  sangerseqr <- requireNamespace("sangerseqR")
  stopifnot(isTRUE(sangerseqr))
  abif <- sangerseqR::read.abif(path)
  if (is.null(abif@data$PCON.2)) {
    message(sprintf("failed on %s", seqname))
    return()
  }
  nucseq <- substring(abif@data$PBAS.2, 1, length(abif@data$PLOC.2))
  if (!typeof(abif@data$PCON.2) == "integer") {
    num_quals <- utf8ToInt(abif@data$PCON.2)[1:length(abif@data$PLOC.2)]
  }
  else {
    num_quals <- abif@data$PCON.2[1:length(abif@data$PLOC.2)]
  }
  if (isTRUE(recall)) {
    recalled <- sangerseqR::makeBaseCalls(sangerseqR::sangerseq(abif))
    nucseq <- sangerseqR::primarySeq(recalled, string = TRUE)
    if (nchar(nucseq) != length(num_quals)) {
      trim <- FALSE
      num_quals <- rep(60, nchar(nucseq))
      warning("Length of quality scores does not equal length of\n              re-called base sequence, ignoring quality scores")
    }
  }
  if (trim == FALSE) {
    tmp1 = list(seqname = seqname, seq = nucseq, 
                quals = rawToChar(as.raw(num_quals + offset)))
    return(tmp1)
  }
  trim_msg <- "Sequence %s can not be trimmed because it is shorter than the trim\n               segment size"
  if (nchar(nucseq) <= min_seq_len) {
    warning(sprintf(trim_msg, seqname))
    return()
  }
  scores = cutoff - 10^(num_quals/-10)
  running_sum <- rep(0, length(scores) + 1)
  for (i in 1:length(scores)) {
    num <- scores[i] + running_sum[i]
    running_sum[i + 1] <- ifelse(num < 0, 0, num)
  }
  trim_start <- min(which(running_sum > 0)) - 1
  trim_finish <- which.max(running_sum) - 2
  if (trim_finish - trim_start < min_seq_len - 1) {
    warning(sprintf(trim_msg, seqname))
    return()
  }
  tmp2 = list(seqname = seqname, seq = substring(nucseq, 
                                                 trim_start, trim_finish), quals = rawToChar(as.raw(num_quals[trim_start:trim_finish] + 
                                                                                                      offset)))
  return(tmp2)
}



trim_ends = function(fastq, threshold){
  
}

### Takes sequence and runs through a while loop until they are properly trimmed
### Appears for CTSS 11 sample
align_and_trim = function(pattern_seq, subject_seq, min_continuity = 15){
  
  gap_length = min_continuity - 1
  raw_pattern = lapply(FUN = rep, X = 1:gap_length, x = "[A-Z]") %>%
    lapply(FUN = paste0, X = ., collapse = "") %>%
    unlist()
  gsub_pattern = raw_pattern %>%
    paste0("^", ., "-|") %>%
    paste0(collapse = "") %>%
    paste0(., raw_pattern %>%
             paste0("-", ., "$|") %>%
             paste0(collapse = ""),
           collapse = "") %>%
    gsub('.{1}$','', .)
  
  input_pattern = pattern_seq
  input_subject = subject_seq
  
  alignment = pairwiseAlignment(pattern = pattern_seq, subject = subject_seq)
  
  output_pattern = alignment@pattern %>% as.character() %>% gsub(gsub_pattern, "", ., perl = T) %>% gsub("-", "", ., perl = T)
  output_subject = alignment@subject %>% as.character() %>% gsub(gsub_pattern, "", ., perl = T) %>% gsub("-", "", ., perl = T)
  
  while(input_pattern != output_pattern | input_subject != output_subject){
    alignment = pairwiseAlignment(pattern = output_pattern, subject = output_subject)
    
    old_output_pattern = output_pattern
    old_output_subject = output_subject
    
    new_output_pattern = alignment@pattern %>% as.character() %>% gsub(gsub_pattern, "", ., perl = T) %>% gsub("-", "", ., perl = T)
    new_output_subject = alignment@subject %>% as.character() %>% gsub(gsub_pattern, "", ., perl = T) %>% gsub("-", "", ., perl = T)
    
    output_pattern = new_output_pattern
    output_subject = new_output_subject
    
    input_pattern = old_output_pattern
    input_subject = old_output_subject
  }
  
  return(list("pattern" = alignment@pattern %>% as.character() %>% gsub("-", "", .),
              "subject" = alignment@subject %>% as.character() %>% gsub("-", "", .),
              "alignment" = alignment))
}

### Substitute multiple string positions simultaneously
subchar <- function(string, pos, char) { 
  for(i in pos) { 
    substr(string, i, i) = char
  } 
  string 
} 

### Make a dataframe with all of the ZAGA information
### Backup ZAGA vector
replacement_zaga = c(rep(0, 989), 0.00998720389310502, 0.00998813447664401,0.009992887520785,
                     0.00999585366068316, 0.00999623914632598, 0.00999799013526835, 0.010001499423723,
                     0.0100030237039207, 0.0100045782875701, 0.0100048452355807, 0.0100049548867042)

make_ZAGA_df = function(sanger_df, p_adjust){
  nvals <- list()
  nvals$A = filter(sanger_df, max_base != "A")$A_area 
  nvals$C = filter(sanger_df, max_base != "C")$C_area 
  nvals$G = filter(sanger_df, max_base != "G")$G_area 
  nvals$T = filter(sanger_df, max_base != "T")$T_area 
  
  n_models =   n_models <-lapply(nvals, FUN = function(x){
    set.seed(1)
    if((unique(x)[1] == 0 & length(unique(x)) == 1) |
       (unique(x)[1] == 0 & length(unique(x)) == 2 & table(x)[2] == 1))
      {x = replacement_zaga; message("Replacement vector used for low noise.")} # add noise if all 0s, or all 0s and one other value.
    tryCatch(gamlss((x)~1, family = ZAGA), error=function(e) # Progressively step up the mu.start if it fails
      tryCatch(gamlss((x)~1, family = ZAGA, mu.start = 1), error=function(e) 
        tryCatch(gamlss((x)~1, family = ZAGA, mu.start = 2), error=function(e) 
          tryCatch(gamlss((x)~1, family = ZAGA, mu.start = 3), error=function(e) # additional step added.
            tryCatch(gamlss((x)~1, family = ZAGA, mu.start = 4), error=function(e) # additional step added.
          gamlss((x)~1, family = ZAGA, mu.start = mean(x))
            )
          )
        )
      )
    )
     # throws errors when a completely 0 vector
  })
  
  null_m_params = lapply(n_models, FUN = function(x){
    mu <- exp(x$mu.coefficients[[1]])
    sigma <- exp(x$sigma.coefficients[[1]])
    nu.logit <- x$nu.coefficients[[1]]
    nu <- exp(nu.logit)/(1+exp(nu.logit))
    fillibens <-cor(as.data.frame(qqnorm(x$residuals, plot = FALSE)))[1,2]
    crit = qZAGA(p = 1-p_adjust, mu = mu, nu = nu, sigma = sigma)
    
    return(data.frame(mu= mu, sigma = sigma, nu = nu, crit = crit, fillibens = fillibens))
  })
  
  null_m_params %>%
    plyr::ldply(., "data.frame") %>%
    dplyr::rename(base = `.id`) %>%
    return()
}

### convert negatives to zeros
neg_to_zero = function(x){ifelse(x < 0, 0, x)}

### Adjust the height and percent values based on GBM and the significance
pvalue_adjust = function(sanger_df, wt, boi, motif, sample_file, critical_values, zaga_parameters, p_value){
  sanger_df %>%
    # only bases of interest
    dplyr::filter(grepl(wt, ctrl_max_base)) %>% # To only pull out the bases of interest from the motif
    dplyr::select(index, ctrl_index, ctrl_max_base, max_base, A_area:T_area, A_perc:T_perc) %>%
    # If the channel is not found in the potential edits, then it is not applicable or NA, otherwise
    # If the height of the channel is greater than the critical value for that channel,
    # return TRUE, otherwise, return FALSE
    mutate(A_pvalue = mapply(FUN = gamlss.dist::dZAGA, x = A_area,
                             mu = zaga_parameters[1, "mu"],
                             sigma = zaga_parameters[1, "sigma"],
                             nu = zaga_parameters[1, "nu"]),
           C_pvalue = mapply(FUN = gamlss.dist::dZAGA, x = C_area,
                             mu = zaga_parameters[2, "mu"],
                             sigma = zaga_parameters[2, "sigma"],
                             nu = zaga_parameters[2, "nu"]),
           G_pvalue = mapply(FUN = gamlss.dist::dZAGA, x = G_area,
                             mu = zaga_parameters[3, "mu"],
                             sigma = zaga_parameters[3, "sigma"],
                             nu = zaga_parameters[3, "nu"]),
           T_pvalue = mapply(FUN = gamlss.dist::dZAGA, x = T_area,
                             mu = zaga_parameters[4, "mu"],
                             sigma = zaga_parameters[4, "sigma"],
                             nu = zaga_parameters[4, "nu"])) %>%
    mutate(A_p_adjust = p.adjust(A_pvalue, "BH"), 
           C_p_adjust = p.adjust(C_pvalue, "BH"), 
           G_p_adjust = p.adjust(G_pvalue, "BH"), 
           T_p_adjust = p.adjust(T_pvalue, "BH")) %>%
    mutate(A_sig = if(!grepl("A", boi)) {FALSE} else {ifelse(A_pvalue <= p_value, TRUE, FALSE)},
           C_sig = if(!grepl("C", boi)) {FALSE} else {ifelse(C_pvalue <= p_value, TRUE, FALSE)},
           G_sig = if(!grepl("G", boi)) {FALSE} else {ifelse(G_pvalue <= p_value, TRUE, FALSE)},
           T_sig = if(!grepl("T", boi)) {FALSE} else {ifelse(T_pvalue <= p_value, TRUE, FALSE)}) %>%
    mutate(A_sig_adjust = if(!grepl("A", boi)) {FALSE} else {ifelse(A_p_adjust<= p_value, TRUE, FALSE)},
           C_sig_adjust = if(!grepl("C", boi)) {FALSE} else {ifelse(C_p_adjust <= p_value, TRUE, FALSE)},
           G_sig_adjust = if(!grepl("G", boi)) {FALSE} else {ifelse(G_p_adjust <= p_value, TRUE, FALSE)},
           T_sig_adjust = if(!grepl("T", boi)) {FALSE} else {ifelse(T_p_adjust <= p_value, TRUE, FALSE)}) %>%
    ### also selecting for the sample index added 07.09.2019
    dplyr::select(index, ctrl_index, ctrl_max_base, max_base, A_area:T_area, A_perc:T_perc, A_sig:T_sig, A_sig_adjust:T_sig_adjust, A_pvalue:T_pvalue, A_p_adjust:T_p_adjust) %>%
    mutate(motif = motif, sample_file = sample_file)
}

### Reverse complements a DNA character string
revcom = function(x){as.character(reverseComplement(DNAString(x)))}

### Plotting functions
plotRawSample = function(raw_sample_df, sample_alt, pre_cross_align_sample_df){
  raw_sample_df %>%
    ggplot(aes(x = index, y = max_base_perc)) +
    scale_y_continuous(limits = c(0,100), breaks = seq(0,100, 10), expand = c(0,0)) +
    scale_x_continuous(limits = c(0, max(raw_sample_df$index)), expand = c(0,0)) +
    geom_bar(data = sample_alt, aes(x = index, y = 100), stat = "identity", fill = "#53BCC2", color = "#53BCC2") +
    geom_rect(xmin = min(pre_cross_align_sample_df$index), ymin = 0,
              xmax = max(pre_cross_align_sample_df$index), ymax = 100,
              fill = "white", color = "black", alpha = 0) +
    xlab("Position in sample file") +
    ylab("Percent signal of basecall") +
    geom_hline(yintercept = mean(pre_cross_align_sample_df$max_base_perc), color = "darkred", size = 1) +
    geom_line() +
    geom_rect(xmin = 0, ymin = 0, xmax = min(pre_cross_align_sample_df$index), ymax = 100,
              fill = "grey", color = "black", alpha = 0.01) +
    geom_rect(xmin = max(pre_cross_align_sample_df$index), ymin = 0,
              xmax = max(raw_sample_df$index), ymax = 100,
              fill = "grey", color = "black", alpha = 0.01) +
    geom_rect(xmin = min(pre_cross_align_sample_df$index) +
                (max(pre_cross_align_sample_df$index) - min(pre_cross_align_sample_df$index))*0.025 - 5,
              ymin = 3,
              xmax = min(pre_cross_align_sample_df$index) +
                (max(pre_cross_align_sample_df$index) - min(pre_cross_align_sample_df$index))*0.025 + 215,
              ymax = 20,
              fill = "white", color = "black") +
    annotate(geom = "text",
             label = paste0("Average percent signal (", round(mean(pre_cross_align_sample_df$max_base_perc), 1),"%)") ,
             color = "darkred",
             x = min(pre_cross_align_sample_df$index) +
               (max(pre_cross_align_sample_df$index) - min(pre_cross_align_sample_df$index))*0.025,
             y = 17,
             hjust = 0,
             size = 6) +
    annotate(geom = "text", label = "Low phred trimmed regions", color = "grey30",
             x = min(pre_cross_align_sample_df$index) +
               (max(pre_cross_align_sample_df$index) - min(pre_cross_align_sample_df$index))*0.025,
             y = 12,
             hjust = 0,
             size = 6) +
    annotate(geom = "text", label = "Motif of interest", color = "#53BCC2",
             x = min(pre_cross_align_sample_df$index) +
               (max(pre_cross_align_sample_df$index) - min(pre_cross_align_sample_df$index))*0.025,
             y = 7,
             hjust = 0,
             size = 6) +
    theme_classic(base_size = 18)
}

# Function for plotting trimmed sample data
plotTrimmedSample = function(sample_df, pre_cross_align_sample_df, output_sample_alt, raw_sample_df, sample_alt){
  sample_df %>%
    mutate(perc = (100-ctrl_max_base_perc) %>% round(., 0)) %>%
    ggplot(aes(x = index, y = perc)) +
    geom_line() +
    scale_y_continuous(limits = c(0,100), breaks = seq(0,100, 10), expand = c(0,0)) +
    scale_x_continuous(limits = c(0, max(raw_sample_df$index)), expand = c(0,0)) +
    # geom_bar(data = sample_alt, aes(x = index, y = 100), stat = "identity", fill = "#53BCC2", color = "#53BCC2", alpha = 0.1) +
    # geom_bar(data = output_sample_alt, aes(x = index, y = 100), stat = "identity", fill = "lightblue", color = "lightblue") +
    geom_rect(xmin = min(pre_cross_align_sample_df$index), ymin = 0,
              xmax = max(pre_cross_align_sample_df$index), ymax = 100,
              fill = "white", color = "black", alpha = 0) +
    xlab("Position in sample file") +
    ylab("Percent noise in WT basecall") +
    geom_hline(yintercept = mean(pre_cross_align_sample_df$max_base_height), color = "darkred", size = 1) +
    geom_line() +
    geom_point(data = output_sample_alt %>%
                 mutate(sig = {ifelse(sig, "Signficant", "Non-significant")}) %>%
                 mutate(perc = (100-ctrl_max_base_perc) %>% round(., 0)),
               aes(fill = motif_id, x = index, y = perc, alpha = sig), pch = 21, size = 2, color = "black") +
    geom_rect(xmin = 0, ymin = 0, xmax = min(pre_cross_align_sample_df$index), ymax = 100,
              fill = "grey", color = "black", alpha = 0.01) +
    geom_rect(xmin = max(pre_cross_align_sample_df$index), ymin = 0,
              xmax = max(raw_sample_df$index), ymax = 100,
              fill = "grey", color = "black", alpha = 0.01) +
    theme_classic(base_size = 18) +
    scale_alpha_manual(values = c(0.3, 1)) +
    theme(legend.position = "none")
}

# Function for producing the table data on editing
calculateEditingData = function(output_sample_alt, edit) {
  output_sample_alt %>%
    dplyr::select(index, ctrl_index, A_perc:T_perc) %>%
    gather(base, perc, A_perc:`T_perc`) %>%
    inner_join(., output_sample_alt %>%
                 dplyr::select(index, A_sig:T_sig) %>%
                 gather(base_sig, sig, A_sig:T_sig)
    ) %>%
    mutate(base = gsub("_perc", "", base), base_sig = gsub("_sig", "", base_sig)) %>%
    mutate(perc = perc*100, sig = {ifelse(sig, "Significant", "Non-significant")}) %>%
    filter(base == base_sig) %>%
    filter(grepl(pattern = edit, x = base)) %>%
    mutate(tmp = 1) %>%
    group_by(base) %>%
    mutate(tally = cumsum(tmp)) %>%
    dplyr::select(-tmp, -base_sig)
}

plotEditingData = function(editing_data) {
  editing_data %>%
    ggplot(aes(x = sig, y = perc, color = sig)) +
    scale_y_continuous(limits = c(0,100), breaks = seq(0,100,10)) +
    geom_bar(aes(fill = sig), stat = "summary", fun.y = "mean",
             color = "black", alpha = 0.3, show.legend = F) +
    geom_jitter(size = 2, alpha = 0.7) +
    ylab("Percent height") +
    xlab("") +
    labs(color = "Potential edit") +
    guides(fill = NULL) +
    theme_classic(base_size = 18) +
    theme(aspect.ratio = 4/1,
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +#sum(grepl(edit, bases))) +
    facet_wrap(.~base, nrow = 1, scales = "free_x")
}

tableEditingData = function(editing_data) {
  editing_data %>%
    group_by(base, sig) %>%
    dplyr::summarize(Mean = mean(perc), Max = max(perc), Min = min(perc), SD = sd(perc), N = length(perc)) %>%
    dplyr::rename(Base = base, Significance = sig)
}

# Define geom_chromatogram function
# input the sanger object, start index, end index and output chromatogram with underlying base percents
geom_chromatogram = function(sanger, start_index, end_index){
  
  # define bases
  bases = c("A", "C", "G", "T")
  
  # define the basecalls  in the sanger object initially to set the proper phasing to the data
  sanger = makeBaseCalls(sanger)
  
  # make rawtraces from the sanger trace matrix
  rawtraces = sanger@traceMatrix %>%
    # the trace matrix is the scatter of points that make up the chromatogram
    # convert the large matrix of traces into a tibble object
    as_tibble() %>%
    # name the columns by base channel
    dplyr::rename(A = V1, C = V2, G = V3, `T` = V4) %>%
    # create a column 1 through n trace measurements taken
    rowid_to_column("x")
  
  # create the peakPosMatrix giving the x coordinate of each basecall
  peakPosMatrix = sanger@peakPosMatrix %>% as_tibble()
  # create column names for each channel
  colnames(peakPosMatrix) = c("A", "C", "G", "T")
  # format peakPosMatrix
  peakPosMatrix = peakPosMatrix %>%
    # add the basecall for each position
    mutate(basecall = sanger@primarySeq %>% as.character %>% strsplit(split = "") %>% unlist()) %>%
    rowid_to_column("index") %>%
    gather(base, x, A:`T`) %>%
    filter(basecall == base) %>%
    arrange(index) %>%
    dplyr::select(index, basecall, x)
  
  # Make the traces data frame
  traces = peakPosMatrix %>% 
    # filter only the trace information at positions that are registered as the position of the peak height
    inner_join(., rawtraces) %>%
    # define the previous base index
    mutate(x_1 = lag(x, default = 1)) %>%
    # start the chromatogram for each base as the position between the two bases divided by two
    # start the chromatogram as the x coordinate bewteen the first  base of interest and the n-1 base
    mutate(start = floor((x + x_1) / 2)) %>%
    # end the chromatogram for each base as the position
    # empirically derived I believe
    mutate(end = lead(start, default = (max(x) + 6))) %>%
    # define the total area for each basecall
    mutate(Tot.Area = A + C + G + `T`) %>%
    # calculate percent heights for each base from the total area
    mutate(A_perc = round(A*100/Tot.Area),
           C_perc = round(C*100/Tot.Area),
           G_perc = round(G*100/Tot.Area),
           T_perc = round(`T`*100/Tot.Area))
  
  # The position list is for every basecall, these are the corresponding x-coordinates
  # This is a list object, where the primary index of the list represents the basecall index
  # The vector element under each index represents the x-coordinates of the peakAmpMatrix for each basecall
  position_list = mapply(FUN = seq,
                         from = c(1, traces$end),
                         to = c(traces$end, max(traces$end)+10),
                         by = 1) %>%
    head(., -1)
  
  # make the names of the position_list indentical to the index of the position list
  names(position_list) = c(1:length(position_list))
  
  indices = ldply(position_list, "data.frame") %>%
    dplyr::rename(index = ".id", x = "data.frame") %>%
    as_tibble() %>%
    mutate(index = as.numeric(index))
  
  # This joins the basecall index to the rawtrace x coordinate based data
  rawtraces %<>% inner_join(., indices)
  
  # Enter the plotting based on the indices
  # defines the basecall_indices as the start and end indices of interest for making the chromatogram
  basecall_indices = c(start_index, end_index)
  
  # determine the start and end x coordinates of the rawtraces to filter on for plotting
  index_filters = rawtraces %>%
    # filter the rawtraces data to just include the indices that are within the region of interest
    # some of the bleedthrough is still observed at this point
    filter(index >= basecall_indices[1] & index <= basecall_indices[2]) %>%
    # the x coordinates corresponding to the region of interest
    .$x %>%
    # takes the min and max of the peakAmpMatrix x cooridinates
    quantile(., c(0, 1)) %>%
    # adds an x coordinate cushion
    {. + c(-6, 6)}
  
  # gather the rawtraces data for ggplot format
  # filter in only traces that fall within the x coordinates rawtraces index filter
  # could merge this line of code with the preceeding index_filters chunk
  plotTraces = rawtraces %>%
    gather(base, height, A:`T`) %>%
    filter(x >= index_filters[1] & x <= index_filters[2]) %>%
    # This join operation is used to creat a base code that later allows a manual gradient to be employed
    inner_join(.,
               data.frame(stringsAsFactors = FALSE,
                          base = c("A", "C", "G", "T"),
                          base_code = c("1000", "2000", "3000", "4000"))
    ) %>%
    # calculate the min and max height for geom_ribbon
    mutate(max_height = pmax(height), min_height = 0)
  
  # calculate the y coordinates for the chromatogram
  y_bases = max(plotTraces$height)*(1/2) %>%
    {.*-c(2/11,4/11,6/11,8/11,10/11)}
  
  # calculate the number of bases involved
  n_bases = basecall_indices[2]-basecall_indices[1]+1
  
  # calculate the x coordinates for the chromatogram
  x_bases = (index_filters[2]-index_filters[1]) %>%
    {.*(seq(2, 2*n_bases+1, 2)/(2*n_bases + 2))} %>%
    {.+index_filters[1]}
  
  tile_plot = traces %>%
    filter(index >= basecall_indices[1] & index <= basecall_indices[2]) %>%
    dplyr::select(basecall, A_perc, C_perc, G_perc, T_perc) %>%
    gather(col, labels, basecall:T_perc) %>%
    mutate(y = rep(y_bases, each = n_bases),
           x = rep(x_bases, times = 5))
  
  base_annotation = data.frame(label = bases,
                               x = rep(min(plotTraces$x), 4),
                               y = tail(rev(sort(unique(tile_plot$y))), -1))
  
  ## Create functions for the fill color palettes
  colors1 = colorRampPalette(colors = c("white", "white"))
  colors2 = colorRampPalette(colors = c("white", "#e41a1c"))
  colors3 = colorRampPalette(colors = c("#e41a1c", "#e41a1c"))
  colors4 = colorRampPalette(colors = c("#e41a1c", "#377eb8"))
  
  # establish a fill_key vector
  fill_key = c(0:6, 7:19, 20:49, 50:100, 1000, 2000, 3000, 4000)
  
  # establish the colors in the fill_key
  fill_colors = c(colors1(length(0:6)),
                  colors2(length(7:19)),
                  colors3(length(20:49)),
                  colors4(length(50:100)),
                  "#32CD32", "#4169E1", "#121212", "#FC4338")
  
  # tie together the fill_key and fill_colors using the names() function
  names(fill_colors) = fill_key
  
  chromatogram = plotTraces %>%
    mutate(max_height = pmax(height), min_height = 0) %>%
    ggplot(data = ., aes(x = x, y = height)) +
    geom_ribbon(aes(ymin = min_height, ymax = max_height, color = base, fill = base_code), alpha = 0.1) +
    scale_color_manual(values = c("A" = "#4daf4a", "C" = "#377eb8", "G" = "black", "T" = "#e41a1c")) +
    geom_tile(data = filter(tile_plot, col != "basecall"),
              aes(y = y, x = x, fill = as.character(as.numeric(labels))), color = "black") +
    scale_fill_manual(values = fill_colors) +
    geom_text(data = tile_plot, aes(x = x, y = y, label = labels), color = "black") +
    geom_text(data = base_annotation, aes(x = x, y = y, label = label), color = "black") +
    theme_void(base_size = 36) +
    theme(legend.position = "none",
          aspect.ratio = 1/1.5)
  
  return(chromatogram)
}

runEditR = function(
    # Set parameters
  sample_name,
  sample_file,
  ctrl_file,
  
  motif = "A", # Use IUPAC notation
  motif_fwd = TRUE,
  
  wt = "A", # Enter wt bases of interest with | separation
  edit = "G", # Enter edit of interest with | separation
  ctrl_seq_fwd = TRUE,
  use_ctrl_seq = TRUE,
  p_value = 0.01,
  phred_cutoff = 0.001, ### 0.0001 seems good.
  trim = TRUE,
  boi = paste0(wt, "|", edit),
  bases = c("A", "C", "G", "T")
  
  # Load sequencing files
  # sample_file = "/Users/kluesner/Desktop/Research/EditR/multiEditR/test_samples/RP272_cdna_wt.ab1"
  # ctrl_file = "/Users/kluesner/Desktop/Research/EditR/multiEditR/test_samples/RP272_cdna_ko.ab1"
  
  #sample_file = "/Users/kluesner/Desktop/Research/EditR/multiEditR/ME1/ME1_CmAG/CmAG_012_RP008_2018\ Sep\ 27.ab1"
  #ctrl_file = "/Users/kluesner/Desktop/Research/EditR/multiEditR/ME1/ME1_CmAG/CmAG/CmAG_001_RP008_2018 Sep 27.ab1"
  
  #sample_file = "/Users/kluesner/Desktop/Research/EditR/multiEditR/test_samples/RP272_cdna_wt.ab1"
  #ctrl_file = "/Users/kluesner/Desktop/Research/EditR/multiEditR/test_samples/RP272_cdna_ko.ab1"
){
  suppressWarnings({
  message("initialized.")
  start = Sys.time()
  # Make sangerseq objects
  # Need to flesh out the TRUE statement branch
  if(use_ctrl_seq)
  { #input_seq = abif_to_fastq(path = ctrl_file, cutoff = phred_cutoff)$seq
    input_seq = read_lines(ctrl_file)[2]
    init_ctrl_seq = input_seq
    # if(ctrl_seq_fwd){} else {init_ctrl_seq = revcom(init_ctrl_seq)}
    ctrl_fastq = list()
    ctrl_fastq$seq = input_seq
    ctrl_df = data.frame(max_base = init_ctrl_seq %>% base::strsplit(., split = "") %>% unlist(),
                         base_call = init_ctrl_seq %>% base::strsplit(., split = "") %>% unlist()) %>%
      mutate(index = 1:NROW(max_base))
  } else 
  {
    # Generate ctrl sanger data frame
    # Generate ctrl primary basecalls
    ctrl_sanger = readsangerseq(ctrl_file)
    ctrl_df = make_ctrl_sanger_df(ctrl_sanger)
    init_ctrl_seq = ctrl_df$base_call %>% paste0(., collapse = "")
    if(ctrl_seq_fwd){} else {init_ctrl_seq = revcom(init_ctrl_seq)}
    ctrl_fastq = abif_to_fastq(path = ctrl_file, cutoff = phred_cutoff)
  }
  
  # IF the sample needs to be reversed, then reverse complement the ctrl_seq
  
  # Make sangerseq object
  # Generate samp sanger data frame
  # Generate samp primary basecalls
  sample_sanger = readsangerseq(sample_file)
  sample_df = make_samp_sanger_df(sample_sanger, init_ctrl_seq)
  init_sample_seq = sample_df$primary_base_call %>% paste0(., collapse = "")
  # Genereate phred scores for ctrl and samp, trimming is built in using mott's algorithm
  sample_fastq = abif_to_fastq(path = sample_file, cutoff = phred_cutoff)
  
  # Align the both the ctrl and samp to their fastq filtered sequences
  # reasonable to assume phred scored sequence will always be smaller than the primary seq
  # use high gap penalty to force SNP alignment
  sample_alignment = pairwiseAlignment(pattern = sample_fastq$seq, subject = init_sample_seq)
  ctrl_alignment = pairwiseAlignment(pattern = ctrl_fastq$seq, subject = init_ctrl_seq)
  
  # Save unfiltered dataframes
  raw_sample_df = sample_df
  raw_ctrl_df = ctrl_df
  
  message("samples loaded.")
  
  # Filter dfs on high phred sequence
  sample_df %<>% 
    filter(index >= sample_alignment@subject@range@start) %<>%
    filter(index <= sample_alignment@subject@range@start + sample_alignment@subject@range@width - 1) %<>%
    mutate(post_filter_index = 1:NROW(index))
  ctrl_df %<>% 
    filter(index >= ctrl_alignment@subject@range@start) %<>%
    filter(index <= ctrl_alignment@subject@range@start + ctrl_alignment@subject@range@width - 1) %<>%
    mutate(post_filter_index = 1:NROW(index))
  
  # Regenerate primary basecalls
  ctrl_seq = ctrl_df$base_call %>% paste0(., collapse = "")
  sample_seq = sample_df$max_base %>% paste0(., collapse = "")
  
  # Save the pre-cross alignment dataframes
  pre_cross_align_sample_df = sample_df
  pre_cross_align_ctrl_df = ctrl_df
  
  message("samples filtered.")
  ### Bring samples together ###
  ### 01.07.19, if a ctrl sequence is used instead of a ctrl sequence, it would enter here.
  ### To use the context correction, would you still need to have the ctrl sequence to apply the GBM to ?
  # Align sample_seq to ctrl_seq
  trimmed_alignment = align_and_trim(sample_seq, ctrl_seq, min_continuity = 15)
  samp_alignment_seq = trimmed_alignment$alignment@pattern %>% as.character()
  ctrl_alignment_seq = trimmed_alignment$alignment@subject %>% as.character()
  
  # Align the trimmed sequences to the sequences from the data frame
  sample_trimmed_alignment = pairwiseAlignment(pattern = trimmed_alignment$pattern, subject = sample_seq, gapOpening = 1000, gapExtension = 1000, type = "local")
  ctrl_trimmed_alignment = pairwiseAlignment(pattern = trimmed_alignment$subject, subject = ctrl_seq, gapOpening = 1000, gapExtension = 1000, type = "local")
  
  ### Add predicted values from GBM
  # sample_df = sample_df %>%
  #   dplyr::select(A_area:T_perc,
  #                 max_base, Tot.Area, index, max_base_height, post_filter_index)
  
  
  
  # Filter dfs to aligned sequences
  sample_df %<>% 
    filter(post_filter_index >= sample_trimmed_alignment@subject@range@start) %<>%
    filter(post_filter_index <= sample_trimmed_alignment@subject@range@start + sample_trimmed_alignment@subject@range@width - 1) %>%
    mutate(post_aligned_index = 1:NROW(index))
  ctrl_df %<>% 
    filter(post_filter_index >= ctrl_trimmed_alignment@subject@range@start) %<>%
    filter(post_filter_index <= ctrl_trimmed_alignment@subject@range@start + ctrl_trimmed_alignment@subject@range@width - 1) %>%
    mutate(post_aligned_index = 1:NROW(index))
  
  message("sample aligned.")
  
  
  ### Generate post-filter, post-aligned sequence
  ctrl_seq = ctrl_df$base_call %>% paste0(., collapse = "")
  samp_seq = sample_df$max_base %>% paste0(., collapse = "")
  
  ### Save a df for NGS analysis
  pre_aligned_sample_df = sample_df
  
  ### Filter out any base positons that have an indel
  samp_indel = samp_alignment_seq %>% gregexpr("-", .) %>% unlist
  ctrl_indel = ctrl_alignment_seq %>% gregexpr("-", .) %>% unlist
  
  ### Code added 10.27.18, as when there is no indels in one sample it returned a -1 instead of 0, which introduced an error
  if(samp_indel == -1){samp_indel = 0} else {}
  if(ctrl_indel == -1){ctrl_indel = 0} else {}
  
  ctrl_df = ctrl_df %>%
    mutate(., 
           indel_filter = ctrl_alignment_seq %>%
             subchar(., samp_indel, "_") %>%
             gsub("-", "", .) %>%
             strsplit(x = ., split = "") %>%
             .[[1]]
    ) %>%
    filter(indel_filter != "_") %>%
    dplyr::select(-indel_filter) %>%
    mutate(filtered_index = 1:NROW(max_base))
  
  sample_df = sample_df %>%
    mutate(., 
           indel_filter = samp_alignment_seq %>%
             subchar(., ctrl_indel, "_") %>%
             gsub("-", "", .) %>%
             strsplit(x = ., split = "") %>%
             .[[1]]
    ) %>%
    filter(indel_filter != "_") %>%
    dplyr::select(-indel_filter) %>%
    mutate(filtered_index = 1:NROW(max_base)) %>% 
    mutate(ctrl_post_aligned_index = ctrl_df$post_aligned_index)
  
  message("Indels removed.")
  ### join the initial ctrl index to the sample to give a reference
  sample_df = ctrl_df %>%
    dplyr::select(post_aligned_index, index) %>%
    dplyr::rename(ctrl_post_aligned_index = post_aligned_index, ctrl_index = index) %>%
    inner_join(., sample_df) #%>%
  #dplyr::select(everything(), ctrl_index, ctrl_post_aligned_index, pre_trinucleotide, post_trinucleotide)
  
  ### Assign the sample and ctrl file names
  sample_df %<>% mutate(sample_file = sample_file, ctrl_file = ctrl_file)
  ctrl_df %<>% mutate(sample_file = sample_file, ctrl_file = ctrl_file)
  
  ### reassign the sequences post filtering
  samp_df_seq = sample_df$max_base %>% paste0(., collapse = "")
  ctrl_df_seq = ctrl_df$max_base %>% paste0(., collapse = "")
  
  message("Indices adjusted.")
  ### ENTER MOTIF ISOLATION ###
  # Will want to be compare this to post alignment seq and index, but before indel removal
  # ctrl_seq will have the same indexing as post_aligned_index, but the indexing is not the same across the two
  # Will need to figure out how to be able to pull the same indexing out of both of them
  # Could take the post_aligned_index from ctrl and apply it to the sample
  
  # Reverse complement motif if needed
  if(motif_fwd){}else{motif = revcom(motif)}
  
  # Align the motif of interest to the ctrl_seq
  motif_alignment = matchPattern(pattern = DNAString(motif), subject = DNAString(ctrl_seq), fixed = FALSE)
  n_alignments = motif_alignment@ranges %>% length()
  
  motif_positions = mapply(FUN = seq,
                           from = motif_alignment@ranges@start,
                           to = (motif_alignment@ranges@start + nchar(motif) - 1)) %>% as.vector()
  names(motif_positions) = rep(x = c(1:n_alignments), each = nchar(motif))
  
  # Append the sequences from the ctrl df to the sample df
  sample_df %<>% mutate(ctrl_max_base = ctrl_df$max_base, ctrl_base_call = ctrl_df$base_call)
  
  # calculate an editing index
  
  # sample_df %<>%
  #   group_by(ctrl_max_base) %<>%
  #   mutate(EI = sum(!! sym(paste0(edit, "_perc"))) / sum(!! sym(edit)) + sum(!! sym(wt)))
  
  # Generate null and alternative samples for distribution
  # Perform for both the sample df
  sample_null = sample_df %>% filter(!(ctrl_post_aligned_index %in% motif_positions)) # This dataframe consists of all bases in the sample where the motif is not found in the ctrl sequence
  sample_alt = sample_df %>% filter(ctrl_post_aligned_index %in% motif_positions) # This dataframe consists of all bases in the sample where the motif is found in the ctrl sequence
  
  # Find all potential events of significant noise
  filtered_sample_alt = sample_alt %>%
    filter(grepl(wt, ctrl_max_base)) %>% # Use the ctrl wt base for determining data
    dplyr::rename(A = A_area, C = C_area, G = G_area, `T` = T_area) %>%
    gather(base, height, A:`T`) %>%
    filter(grepl(edit, base)) # Filter out hypothetical mutations that are not of interest
  
  message("Motifs of interest mapped and subsetted.")
  
  # Adjust p-value
  n_comparisons = NROW(filtered_sample_alt)
  # Holm-sidak correction, may need to use the smirnov correction for family-wise error rates
  # Will need to read original paper and cite appropriately
  # if(adjust_p) {p_adjust = p.adjust(p_value, method = "holm", n = n_comparisons)} else {p_adjust = p_value}
  
  # Generate zG models for each base
  # uses the sample_null to calculate
  zaga_parameters = make_ZAGA_df(sample_null, p_adjust = p_value) %>%
    mutate(sample_file = sample_file) %>%
    mutate(sample_name = sample_name) %>%
    dplyr::select(sample_name, everything())
  
  critical_values = zaga_parameters$crit
  
  ### Find significant edits and then apply the GBM adjustments
  # Determine which values are significant
  # Keep significant values and replace all n.s. values with NA
  # Use the significant values to calculated an adjusted height using formula 1
  # Find all significant edits
  output_sample_alt = pvalue_adjust(sample_alt, wt, boi, motif, sample_file, critical_values, zaga_parameters, p_value)

  
  # Create Editing index for output_sample_alt
  # 02.02.2020
  output_sample_alt = output_sample_alt %>%
    dplyr::group_by(ctrl_max_base) %>% # technically shouldn't change anything as it's the same across samples
    dplyr::mutate(AEI_sanger = (sum(G_perc) / (sum(A_perc) + sum(G_perc))),
                  CEI_sanger = (sum(T_perc) / (sum(C_perc) + sum(T_perc))),
                  GEI_sanger = (sum(A_perc) / (sum(G_perc) + sum(A_perc))),
                  TEI_sanger = (sum(C_perc) / (sum(T_perc) + sum(C_perc)))
    ) %>%
    ungroup() %>%
    # Keep only the EI_index for each reference base
    gather(EI_base, EI_sanger, AEI_sanger:TEI_sanger) %>%
    mutate(EI_base = gsub("EI_sanger", "", EI_base)) %>%
    filter(EI_base == ctrl_max_base) %>%
    dplyr::select(-EI_base)
  
  # define control chromatogram indices to set base position
  sample_chromatogram_indices = range(output_sample_alt$index) %>% sort
  ctrl_chromatogram_indices = range(output_sample_alt$ctrl_index) %>% sort
  
  output_sample = output_sample_alt %>%
    mutate(target_base = ctrl_index - ctrl_chromatogram_indices[1] + 1,
           sample_name = sample_name) %>%
    dplyr::select(sample_name, target_base, `motif`, ctrl_max_base, A_perc:T_perc, A_sig:T_sig, A_pvalue:T_pvalue, index, ctrl_index, sample_file)
  
  # output_sample_data = output_sample_alt %>%
  #   mutate(base_position = ctrl_index - ctrl_chromatogram_indices[1] + 1) %>%
  #   dplyr::select(base_position, motif, index, ctrl_index, ctrl_max_base, A_perc:T_perc, A_sig:T_sig, A_pvalue:T_pvalue, sample_file)
  # 
  # 
  output_sample_null = sample_null %>%
    dplyr::select(ctrl_index, ctrl_max_base, max_base, A_area:T_area, A_perc:T_perc) %>%
    mutate(motif = motif, sample_file = sample_file)
  
  ### Plot chromatograms
  sample_chromatogram = geom_chromatogram(sample_sanger, start = sample_chromatogram_indices[1], end = sample_chromatogram_indices[1] + nchar(motif)-1) +
    ggtitle(paste0(sample_name, " sample")) +
    theme(plot.title = element_text(size = 16))
  
  
  if(use_ctrl_seq){
    ctrl_chromatogram = NULL
  } else
  {
    ctrl_chromatogram = geom_chromatogram(ctrl_sanger, start = ctrl_chromatogram_indices[1], end = ctrl_chromatogram_indices[1] + nchar(motif)-1) +
      ggtitle(paste0(sample_name, " control")) +
      theme(plot.title = element_text(size = 16))
  }
  
  # chromatograms = grid.arrange(arrangeGrob(sample_chromatogram, ctrl_chromatogram, ncol = 1, nrow = 2), heights = 20, widths = 30)

  
  ### 
  output = list(#"sample_data" = sample_df,
                "sample_data" = output_sample,
                # "control_sequence" = init_ctrl_seq,
                # "control_data" = ctrl_df, 
                "statistical_parameters" = zaga_parameters,
                "sample_chromatogram" = sample_chromatogram,
                "control_chromatogram" = ctrl_chromatogram
                )
  
  # Reset start directory
  message(paste0(round(Sys.time() - start, 2), " seconds elapsed."))
  
  return(output)
  })
}

# runEditR in batch mode
runEditRBatch = function(row, params){
  message(paste0(params[row,2] %>% unlist, " started."))
  runEditR(sample_name = params[row,2] %>% unlist,
           sample_file = params[row,3] %>% unlist,
           ctrl_file = params[row,4] %>% unlist,
           motif = params[row,5] %>% unlist,
           motif_fwd = params[row,6] %>% unlist,
           wt = params[row,7] %>% unlist, 
           edit = params[row,8] %>% unlist,
           ctrl_seq_fwd = params[row,9] %>% unlist,
           use_ctrl_seq = params[row,10] %>% unlist,
           p_value = params[row,11] %>% unlist)
}
