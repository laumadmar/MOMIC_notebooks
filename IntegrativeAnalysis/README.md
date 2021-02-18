# Integrative Analysis tutorial

## Short manual on how to perform your own integrative analysis. Read user manual to get more details

Integrative Analysis aims at combining heterogeneous data at different omic levels. 

The template is provided as read-only. To create your own notebook from the template, select it and click on the Duplicate button. Modify paths and code according to your needs. Alternatively, create an empy notebook from the New button located at the top left menu and folow the same steps indicated in the template.

Execute code in cells with `CTRL+ENTER` or doing click on the Run button located at the top menu. Alternative `CTRL+ALT` to create a new empty code cell bellow it. Click on Help button-> Keyboard Shortcuts to know some useful shortcuts. 

The integration is performed using Robust Rank Aggregation method. It takes as input two or more lists of genes ordered by pvalue ascending (or whichever criteria you consider to order the input genes) and detects genes that are ranked consistently better than expected under null hypothesis of uncorrelated inputs and assigns a significance score for each gene. Final gene list is ranked according to this global score, using the Rank library from the R Basic package. The same rank is assign to those genes with NA score.