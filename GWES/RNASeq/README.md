# RNASeq tutorial 

## Short tutorial on how to perform your own RNAseq Analysis. Read user manual to get more details

RNA-Seq is a particular technology-based sequencing technique which uses next-generation sequencing (NGS) to reveal the presence and quantity of RNA in a biological sample at a given moment.

Notebooks are provided as read-only templates and are self-explanatory. To create your own notebook from a template, select the desired notebook/template and click on the Duplicate button. Modify paths and code according to your needs. Alternatively, create an empty notebook from the New button located at the top left menu and follow the same steps indicated in the template.

Execute code in cells with `CTRL+ENTER` or doing click on the Run button located at the top menu. Alternative `CTRL+ALT` to create a new empty code cell bellow it. Click on Help button-> Keyboard Shortcuts to know some useful shortcuts. 

A typical workflow to perform RNASeq analysis is as follow:
1. QC. Check the overall quality of the sequenced reads with fastqc tool usign the template Task1_QC_raw_data
2. Alignment. Align the reads to the reference genome with STAR using the template Task2.1_Alignment_and_ReadQuantification. Run Task2.2_ReadQuantification_to_DESeqDataSet to transform the results of STAR quantification into a format that can be used in the DE step.
3. Differential Expression. Perform DE analysis with DESeq2 using the template Task3_DifferentialAnalysis_compact or the step_by_step one for a detailed protocol without functions wrapping the code.

