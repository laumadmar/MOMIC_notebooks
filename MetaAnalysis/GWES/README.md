# GWES Meta Analysis tutorial

## Short manual on how to perform your own meta analysis. Read user manual to get more details

This meta analysis combines the results from differential expression analysis of distinct datasets, generating average effect sizes and p values for one gene across the different datasets. Genes will be ranked according to this global p value, using the Rank library from the R Basic package. 

The template is provided as read-only. To create your own notebook from the template, select it and click on the Duplicate button. Modify paths and code according to your needs. Alternatively, create an empy notebook from the New button located at the top left menu and follow the same steps indicated in the template.

Note the cell that imports the core pipeline functions `source("/home/guess/scripts/metaDE.R")`. It is very important that you run this cell first.

Execute code in cells with `CTRL+ENTER` or doing click on the Run button located at the top menu. Alternative `CTRL+ALT` to create a new empty code cell bellow it. Click on Help button-> Keyboard Shortcuts to know some useful shortcuts. 

The MetaDE R library (slightly modified) is used to perform this analysis for gene expression data. The algorithm takes as input the p values, observed effect size (logFC values) and observed variance or just p value if the meta analysis is to be done just combining p values. The MetaAnalysis_template.ipynb on this directory combines the two limma tables obtained from the Microarray DE analysis performed for GSE48350 and GSE15222 following the templates on the GWES/Microarray directory.  The core pipeline function used in this template, prepare_matrix_function, takes as input limma datasets containing logFC, confidence intervals and p value. Any output from a different DE algorithm can be used as long as you prepare the datasets accordingly to run the MetaDE algorithm.

Follow the steps on the template and read it carefully to know more about functions and parameters that can be used to customise your analysis. It is also possible to modify a core pipeline function, just paste the function source code in your template, modify it and execute it with the same name. An alternative is to name it differently and use this one instead. 