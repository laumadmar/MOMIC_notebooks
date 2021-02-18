# GWAS Meta Analysis tutorial

## Short manual on how to perform your own integrative analysis. Read user manual to get more details

Different GWAS results can be combined by means of the standard error based weights meta-analysis method implemented in METAL, corrected by study-specific inflation factors. 

The template is provided as read-only. To create your own notebook from the template, select it and click on the Duplicate button. Modify paths and code according to your needs. Alternatively, create an empy notebook from the New button located at the top left menu and folow the same steps indicated in the template.

METAL is a command line tool, so the first cell to be run is `%load_ext rpy2.ipython`, in order to be able to execute bash script.

Execute code in cells with `CTRL+ENTER` or doing click on the Run button located at the top menu. Alternative `CTRL+ALT` to create a new empty code cell bellow it. Click on Help button-> Keyboard Shortcuts to know some useful shortcuts. 