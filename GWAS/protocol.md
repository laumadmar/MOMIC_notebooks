# Complete protolol for GWAS data

An initial standard QC needs to be performed, it consists on: 

    exclusion of individuals: a) with more than 3% missing genotypes, b) with excess autosomal heterozygosity (>0.35 or more than  ± 3 standard deviations from the mean), c) showing a discrepancy between genotypic and reported sex or d) showing non-European ancestry. Duplicated and related individuals are also removed by means of IBS estimates within and across studies. 

    exclusion of SNPs: a) with missing genotype rate > 5%, b) with significant missingness test between cases and controls (p<10−6), c) not in Hardy-Weinberg equilibrium (p<10-6 in controls) and d) SNPs with minor allele frequency (MAF) < 1%. 

After the QC, genotype imputation can be performed with the minimac 3 algorithm at the University of Michigan server using the HRC reference panel, and the SHAPEIT tool for haplotype phasing. After imputation, SNPs with an R2 quality estimate lower than 0.3 is excluded from further analyses according to the software recommendations. 

PLINK can then be used to perform a case control association study.

If you have different datasets to be combined, it can be done by means of the standard error based weights meta-analysis method implemented in METAL, corrected by study-specific inflation factors. A Cochran’s Q test for heterogeneity and I2 estimates is generated to evaluate the potential effect of study heterogeneity on the results.  

Gene-wise statistics can be computed using MAGMA software, which takes into account physical distance and linkage disequilibrium (LD) between markers.

Last, results are inspected usign Manhattanand QQ plots. 

Protocols followed in this pipeline are:
- Anderson, C. A., Pettersson, F. H., Clarke, G. M., Cardon, L. R., Morris, A. P., & Zondervan, K. T. (2010). Data quality control in genetic case‐control association studies
- Marees AT, de Kluiver H, Stringer S, et al. A tutorial on conducting genome-wide association studies. 

The advantage is that all steps are kept in the same notebook, including commands that are traditionally executed in a terminal, so the user can always see past outputs and exact parameters. It helps reproducibility.