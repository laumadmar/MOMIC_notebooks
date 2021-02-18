## Short tutorial on how to perform your own analysis. Read user manual to get more details

GWAS is an observational study of a genome-wide set of genetic variants in different individuals to examine if any variant is associated with a trait.

There are several templates for a complete GWAS analysis using the publicly available 1000 Genome dataset. These are written mainly in bash using the python kernel.

Templates are provided as read-only and are self-explanatory. To create your own notebook from a template, select the desired template and click on the Duplicate button. Modify paths and code according to your needs. Alternatively, create an empty notebook from the New button located at the top left menu and follow the same steps indicated in the template.

Execute code in cells with CTRL+ENTER or doing click on the Run button located at the top menu. Alternative CTRL+ALT to create a new empty code cell bellow it. Click on Help button-> Keyboard Shortcuts to know some useful shortcuts.

Templates provided are for:
- Task1_preQC_template.ipynb: for preparing a working dataset in PLINK v1.9 binary format with all SNPs identified by the rs number and coordinates based on the genome build GRCh37/hg19 (as required by the Michigan imputation server). 
- Task2_QC_template.ipynb: implements a quality control process aimed at removing individuals and markers with particularly high error rates. This QC will be first conducted on a “per-individual” basis prior to QC on a “per-marker” basis, following the guidelines in Anderson et al. protocol.
- Task3_Imputation_template.ipynb: implements genotype imputation. It is performed with the minimac 3 algorithm at the University of Michigan server using the HRC reference panel, and the SHAPEIT tool for haplotype phasing. After imputation, SNPs with an R2 quality estimate lower than 0.3 is excluded from further analyses according to the software recommendations. 
- Task4_Assoc_template.ipynb: performs a case control association study using PLINK. Includes MAGMA.
- Task5_Visualisation_template.ipynb: displays GWAS results, plotting p-values that indicate the significance of the difference in frequency of the allele tested between cases and controls.  
- Task6_Gene-wise_Statistics_template: performs gene-wise statistics woth MAGMA.
