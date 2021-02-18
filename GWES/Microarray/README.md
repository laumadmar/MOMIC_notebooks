# Microarray tutorial

## Short tutorial on how to perform your own Microarray Analysis of Gene Expression. Read user manual to get more details

DNA microarrays can simultaneously measure the expression level of thousands of genes within a particular mRNA sample, which can be used to compare the level of gene transcription in clinical conditions in order to find differences in expression levels between predefined groups of samples.

Notebooks are provided as read-only templates and are self-explanatory. To create your own notebook from a template, select the desired notebook/template and click on the Duplicate button. Modify paths and code according to your needs. Alternatively, create an empty notebook from the New button located at the top left menu and follow the same steps indicated in the template.

Note the cell that imports the core pipeline functions, `source("/mnt/data/GeneExpression/scripts/diffExpressionPipeline.R")`, it is very important that you run this cell first.

Execute code in cells with `CTRL+ENTER` or doing click on the Run button located at the top menu. Alternative `CTRL+ALT` to create a new empty code cell bellow it. Click on Help button-> Keyboard Shortcuts to know some useful shortcuts. 

Before running a Differential Expression analysis, the data (expression and clinical/metadata datasets) have to be clean, containing the same samples and in the same order and ready to use in limma; expression values have to be in logarithm scale. Do a previous pre-processing and raw to intensity matrix transformation if needed. Templates for these tasks are also provided.

Two notebook templates are provided, called compact and step_by_step. Both implement the same functionality but the last one does it in a more detailed and comprehensive way but more tedious to programmatically execute the analysis.

Follow the steps on the template and read it carefully to know more about functions and parameters that can be used to customised your analysis. It is also possible to modify a core pipeline function, just paste the function source code in your template, modify it and execute it with the same name. An alternative is to name it differently and use this one instead.