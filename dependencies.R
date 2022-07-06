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
install.packages("rmarkdown")
install.packages("plotly")
install.packages("yaml")



# Bioconductor packages
# Updated 4.7.19 due to error with Bioconductor packages
# Updated again 4.30.19

BiocManager::install("Biostrings")
BiocManager::install("sangerseqR")
