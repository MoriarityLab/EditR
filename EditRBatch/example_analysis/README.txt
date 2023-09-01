README.txt

Welcome! This Rmarkdown file allows for the batch analysis of base editing from sanger sequencing files

  - To use set up a directory structured like the `example_analysis` directory
      - `example_analysis`
          - `README.txt`
          - `multiEditRBatch.Rmd`
          - `global.R`
          - `parameters.xlsx`
          - `files`
          - `.ab1` or `.fasta` files to be used in the analysis
          - `output`
              - `example_analysis_analysis_data.tsv`, `example_analysis_model_fit.tsv`, and `analysis.html` will be created after running `EditRBatch.Rmd`
          
  - Do not change the directory structure, as it will disrupt the analysis!
      
How to use:

1) Make a copy of the parent `example_analysis` directory, and rename as relevant
     - for e.g. rename to `MGK001`, etc.
     
2) Open `global.R` and install all the packages under libraries if not already installed

3) Open `parameters.xlsx`, edit and save the input data
    - `sample_name` column will be used to identify data in analysis
    - For file names, just use the terminal file name as they'd appear in the final directory
    
4) Move all files referenced in `parameters.xlsx` to the `files` directory

5) Open `EditRBatch.Rmd` in RStudio

6) Edit `analysis_directory = "SYSTEM_PATH"` to reflect the system path to the parent analysis directory

7) Knit `EditRBatch.Rmd`

8) Open `example_analysis/output/analysis.html` to observe results

Disclaimer: This is a preliminary script, with virtually no error handling, so extensive trouble shooting may be required
	
	
