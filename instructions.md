# EditR
 
The goal of this software is to predict where a base edit occurred. Simply, provide your guide RNA protospacer sequence (~ 20bp) and a .ab1 Sanger sequencing file of the edited region (~300-700bp). EditR will then generate distributions of the expected level of noise based on your sample. From these expectations, EditR assigns a probability that each background peak in your guide region is from noise -- in other a words a *P*-value for it being a product of base editing as opposed to noise.

## How to use EditR

1. Upload your .ab1 file of the edited region
2. Enter your gRNA sequence
	+ If your gRNA is antisense to the .ab1 file, check the "Guide sequence is reverse complement" box!
4. On the top of the page, click the "Data QC" tab
5. Check the "Data QA: Signal and noise plot" to see if your sample was trimmed to a region of relatively homogenous noise (orange peaks)
	+ If the noise looks relatively heterogenous, cursor over the plot to see the indices where the homogenous region starts and stops (e.g. 56 and 405)
	+ Enter these values into the left pane under "5' start" and "3' end"
6. If the noise is relatively homogenous, on the top of the page, click the "Predicted Editing" tab
7. Examine the gRNA protospacer chromatogram and underlying tile plot to determine if base editing occured
	+ All colored tiles represent base calls that are deemed significant, i.e. if there multiple colored tiles under a single basecall, base editing likely occurred.
8. If you wish to download a report of the operations performed on your data, on the top of the page click the "Download Report" tab.

You can also provide a custom *P*-value cutoff for calling base-editing, as well as custom ranges to filter your Sanger data on. Check out the QA tab, and scroll down to the plot of post-filtering signal / noise to see if you need to specify custom ranges to filter your data on.

### Output

**Data QC tab**: This tab shows you the quality of the Sanger data that you provided. First is a plot of the total peak area before filtering. Usually the beginning and the end are low quality, so you might want to filter that off. We have a basic filtering algorithm where we cut of the first 20 bases, and then remove any bases that have a total peak area less than a tenth of the mean. The next plot is an interactive plot of the signal and the noise. You can mouse over the plot to get values and the position of the Sanger sequence to use for specifying a custom trimming region. The next plot is a representation of signal / noise as a percent, and the final one is a chromatograph of the specified guide region.

**Predicted Editing tab**: The first plot is a plot of the the peak areas of the different bases along the guide sequence. Boxes are colored in if they are the main peak, or are less than the p-value cutoff for base editing. The second plot is a group of four plots, one for each base, showing the levels of noise at each position, the cutoff for calling a background peak drawn as a line. If the line is much higher than 10, your Sanger sequencing file is likely noisy -- you can check this out on the **Data QC** tab

**Download Report tab:** This tab has a report to download. It contains the predicted editing output, as well as the Data QC plots. It might take some time to compile.

### Special thanks to / sources

Hill JT, Demarest BL, Bisgrove BW, Su YC, Smith M, Yost HJ. (2014) **Poly Peak Parser: Method and software for identification of unknown indels using Sanger Sequencing of PCR products.** *Developmental Dynamics.* [PMID: 25160973](http://www.ncbi.nlm.nih.gov/pubmed/25160973)

Brinkman EK, Chen T, Amendola M, van Steensel B. (2014) **Easy quantitative assessment of genome editing by sequence trace decomposition.** *Nucleic Acids Research.* doi: 10.1093/nar/gku936

## Citation

If our software helps you -- please cite us!

Kluesner MG, Nedveck DA, Lahr WS, Garbe J, Abrahante J, Webber B, Moriarty BS.  **EditR: A method to quantify base editing via Sanger sequencing**. *The CRISPR Journal*. *2018*.

## Troubleshooting

* gRNA sequence is not matching
	+ Firstly, try reverse complement matching the gRNA sequence by checking the "Guide sequence is reverse complement" on the left pane
	+ The .ab1 file may be overly noisy, check the "Data QC" tab first, if so you may need to re-sequence your sample
	+ If there are SNPs within your gRNA region, noisy sequencing or other artifacts, replace ambiguous bases with "N"s in your gRNA sequence
		- EditR can use any IUPAC nucleotides for the gRNA input
		- Additionally, the gRNA sequence length can be extended indefinitely to subset a continous region of interest (e.g. a 30-50 bp region)
* When uploading your .ab1 file the error `<!DOCTYPE HTML><html_lang="en"><head><title>Error</title><style type="text/css">` is given
	+ Open the .ab1 file in a .ab1 file viewer such as SnapGene Viewer
	+ Manually delete the 5' and 3' end of the file to create a region of relatively clean sequencing on either side
	+ Save this edited file as a new .ab1 file, then re-upload to EditR
* Your sequencing is too noisy to detect an edit
	+ Try gel extracting your PCR product
		- Before sending your PCR product for sequencing load 25-50 uL of PCR reaction on a 1% agarose gel in TAE buffer
		- Run the gel and excise your band, and then purify using a gel extraction kit (Qiagen QIAquick works well)
		- Elute your sample in water and send to a trusted Sanger sequencing provider (ACGTinc, single pass works well)
* Your sequencing is too clean to detect an edit
	+ While this may seem paradoxical, EditR actually needs some level of noise for the algorithm to work
	+ To introduce noise that was automotically trimmed out, click on the "Data QC" tab at the top of the page.
	+ See the length of the raw .ab1 file and enter the values of the full length into the "5' start" and "3' end" inputs on the left pane
	+ E.g. enter the values 1 and 452
