# Enrichment analysis tutorial

## Short tutorial on how to perform your own enrichment analysis. Read user manual to get more details

Enrichment analysis, or pathway analysis, can identify terms which are statistically over or under-represented within the list of interest.

The enrichment is performed using WebGestalt in R and the visualization is done in two different representation ways with GOplot package and pheatmap from R CRAN. To ilustrate this we use the results of a meta analysis completed following the GWES pipeline.

Templates are provided as read-only. To create your own notebook from a template, select the desired template and click on the Duplicate button. Modify paths and code according to your needs. Alternatively, create an empty notebook from the New button located at the top left menu and follow the same steps indicated in the template.

Execute code in cells with `CTRL+ENTER` or doing click on the Run button located at the top menu. Alternative `CTRL+ALT` to create a new empty code cell bellow it. Click on Help button-> Keyboard Shortcuts to know some useful shortcuts. 

Templates provided are:
-enrichment_GOplots_template.ipynb: performs the enrichment for all predefined categories and prints the GO plots. 
-enrichment_GOplots_loop_template.ipynb: same as the previous one but includes a loop through more than one dataframe to perform the enrichment on.
-pheatmap_template.ipynb: uses the R library pheatmap to represent combined results from meta expression, meta GWAS integrative analysis and enrichment of meta expression, completed following the different analysis available in this pipeline.


