# run these lines one at a time, and be careful to note if there are 
# any errors put out during installation

# CRAN packages
install.packages("shiny")
install.packages("gamlss")
install.packages("magrittr")
install.packages("dplyr")
install.packages("tidyr")
install.packages("ggplot2")
install.packages("cowplot")
install.packages("gridExtra")
install.packages("rmarkdown") #07.06.2022 if errors that there was no markdown package, try 'install.packages("markdown")'
install.packages("plotly")
install.packages("yaml")



# Bioconductor packages
# Updated 4.7.19 due to error with Bioconductor packages
# Updated again 4.30.19

# Run this code chunk in terminal first
#
# options(repos = BiocManager::repositories())
# source("https://bioconductor.org/biocLite.R")
#
# rsconnect::deployApp("/Users/kluesner/Desktop/Research/spliceR/apps/SpliceR")

BiocManager::install("Biostrings")
BiocManager::install("sangerseqR")
