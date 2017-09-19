## What this is

The goal of this software is to predict where base editing is occurring in a specified region. You provide the guide RNA (usually ~ 20bp), and the output of a Sanger sequencing run of your region of interest, and the prediction software generates an expectation of the noise inherent in your sample. It then takes these expectations, and assigns a probability that each background peak in your guide region is from noise -- essentially a p-value for it being a product of base editing.

### Citation

If our software helps you -- please cite us!

Kluesner MG, Nedveck DA, Lahr WS, Moriarty BS.  **Simple and cost effective quantitative assessment of base editing by quantification of Sanger trace fluorescence.** *Nucleic Acids Research.* *In review*

## Running the prediction software

To run this software, please provide:

*   A sanger sequencing file of your region of interest
*   The guide RNA sequence used for targeted base editing

You can also provide a custom P-value cutoff for calling base-editing, as well as custom ranges to filter your Sanger data on. Check out the QA tab, and scroll down to the plot of post-filtering signal / noise to see if you need to specify custom ranges to filter your data on.

If you want to see some output of the software, just check the “Load example data” checkbox.

### Output

**Data QC tab**: This tab shows you the quality of the Sanger data that you provided. First is a plot of the total peak area before filtering. Usually the beginning and the end are low quality, so you might want to filter that off. We have a basic filtering algorithm where we cut of the first 20 bases, and then remove any bases that have a total peak area less than a tenth of the mean. The next plot is an interactive plot of the signal and the noise. You can mouse over the plot to get values and the position of the Sanger sequence to use for specifying a custom trimming region. The next plot is a representation of signal / noise as a percent, and the final one is a chromatograph of the specified guide region.

**Predicted Editing tab**: The first plot is a plot of the the peak areas of the different bases along the guide sequence. Boxes are colored in if they are the main peak, or are less than the p-value cutoff for base editing. The second plot is a group of four plots, one for each base, showing the levels of noise at each position, the cutoff for calling a background peak drawn as a line. If the line is much higher than 10, your Sanger sequencing file is likely noisy -- you can check this out on the **Data QC** tab

**Download Report tab:** This tab has a report to download. It contains the predicted editing output, as well as the Data QC plots. It might take some time to compile.

### Special thanks to / sources

Hill JT, Demarest BL, Bisgrove BW, Su YC, Smith M, Yost HJ. (2014) **Poly Peak Parser: Method and software for identification of unknown indels using Sanger Sequencing of PCR products.** *Developmental Dynamics.* [PMID: 25160973](http://www.ncbi.nlm.nih.gov/pubmed/25160973)

Brinkman EK, Chen T, Amendola M, van Steensel B. (2014) **Easy quantitative assessment of genome editing by sequence trace decomposition.** *Nucleic Acids Research.* doi: 10.1093/nar/gku936

