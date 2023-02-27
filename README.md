README
================
Yue Pan
27 February, 2023

# Overview

mcQTL (Multi-omic and Cell-type-specific Quantitative Trait Loci) is a
tool that estimates cell type proportions in bulk proteomes by either
using single data source reference or borrowing information in matched
transcriptomes. Based on the deconvoluted cellular composition, mcQTL
further performs Quantitative Trait Loci mapping at cellular resolution,
as well as integrates and visualizes multi-source profiles at bulk and
cell type levels.

# Installation

mcQTL depends on the following packages:

r Biocpkg(“SummarizedExperiment”), for data manipulation,

r Biocpkg(“BiocParallel”), for parallel computing implementation,

r Biocpkg(“TOAST”), for parallel computing implementation.

Please install from github repository “YuePan027/mcQTL” under “dev”
branch.

``` r
devtools::install_github("YuePan027/mcQTL@dev") # from dev branch
#> utf8        (1.2.2  -> 1.2.3 ) [CRAN]
#> fansi       (1.0.3  -> 1.0.4 ) [CRAN]
#> stringi     (1.7.8  -> 1.7.12) [CRAN]
#> cpp11       (0.4.2  -> 0.4.3 ) [CRAN]
#> tidyselect  (1.1.2  -> 1.2.0 ) [CRAN]
#> stringr     (1.4.1  -> 1.5.0 ) [CRAN]
#> dplyr       (1.0.10 -> 1.1.0 ) [CRAN]
#> sys         (3.4    -> 3.4.1 ) [CRAN]
#> openssl     (2.0.2  -> 2.0.5 ) [CRAN]
#> jsonlite    (1.8.0  -> 1.8.4 ) [CRAN]
#> curl        (4.3.2  -> 5.0.0 ) [CRAN]
#> httr        (1.4.4  -> 1.4.5 ) [CRAN]
#> ps          (1.7.1  -> 1.7.2 ) [CRAN]
#> fastmap     (1.1.0  -> 1.1.1 ) [CRAN]
#> digest      (0.6.29 -> 0.6.31) [CRAN]
#> cachem      (1.0.6  -> 1.0.7 ) [CRAN]
#> sass        (0.4.2  -> 0.4.5 ) [CRAN]
#> tinytex     (0.41   -> 0.44  ) [CRAN]
#> bslib       (0.4.0  -> 0.4.2 ) [CRAN]
#> xfun        (0.32   -> 0.37  ) [CRAN]
#> yaml        (2.3.5  -> 2.3.7 ) [CRAN]
#> highr       (0.9    -> 0.10  ) [CRAN]
#> evaluate    (0.16   -> 0.20  ) [CRAN]
#> processx    (3.7.0  -> 3.8.0 ) [CRAN]
#> rmarkdown   (2.16   -> 2.20  ) [CRAN]
#> knitr       (1.40   -> 1.42  ) [CRAN]
#> fs          (1.5.2  -> 1.6.1 ) [CRAN]
#> callr       (3.7.2  -> 3.7.3 ) [CRAN]
#> crayon      (1.5.1  -> 1.5.2 ) [CRAN]
#> bit         (4.0.4  -> 4.0.5 ) [CRAN]
#> vroom       (1.5.7  -> 1.6.1 ) [CRAN]
#> tidyr       (1.2.0  -> 1.3.0 ) [CRAN]
#> broom       (1.0.1  -> 1.0.3 ) [CRAN]
#> timechange  (0.1.1  -> 0.2.0 ) [CRAN]
#> readr       (2.1.2  -> 2.1.4 ) [CRAN]
#> forcats     (0.5.2  -> 1.0.0 ) [CRAN]
#> gargle      (1.2.1  -> 1.3.0 ) [CRAN]
#> colorspace  (2.0-3  -> 2.1-0 ) [CRAN]
#> isoband     (0.2.5  -> 0.2.7 ) [CRAN]
#> data.table  (1.14.2 -> 1.14.8) [CRAN]
#> pkgload     (1.3.0  -> 1.3.2 ) [CRAN]
#> testthat    (3.1.4  -> 3.1.6 ) [CRAN]
#> gtools      (3.9.3  -> 3.9.4 ) [CRAN]
#> formatR     (1.12   -> 1.14  ) [CRAN]
#> readxl      (1.4.1  -> 1.4.2 ) [CRAN]
#> ragg        (1.2.2  -> 1.2.5 ) [CRAN]
#> modelr      (0.1.9  -> 0.1.10) [CRAN]
#> lubridate   (1.8.0  -> 1.9.2 ) [CRAN]
#> ggplot2     (3.3.6  -> 3.4.1 ) [CRAN]
#> dtplyr      (1.2.2  -> 1.3.0 ) [CRAN]
#> dbplyr      (2.2.1  -> 2.3.1 ) [CRAN]
#> pbapply     (1.5-0  -> 1.7-0 ) [CRAN]
#> matrixStats (0.62.0 -> 0.63.0) [CRAN]
#> tidyverse   (1.3.2  -> 2.0.0 ) [CRAN]
#> e1071       (1.7-11 -> 1.7-13) [CRAN]
#> package 'utf8' successfully unpacked and MD5 sums checked
#> package 'fansi' successfully unpacked and MD5 sums checked
#> package 'stringi' successfully unpacked and MD5 sums checked
#> package 'cpp11' successfully unpacked and MD5 sums checked
#> package 'tidyselect' successfully unpacked and MD5 sums checked
#> package 'stringr' successfully unpacked and MD5 sums checked
#> package 'dplyr' successfully unpacked and MD5 sums checked
#> package 'sys' successfully unpacked and MD5 sums checked
#> package 'openssl' successfully unpacked and MD5 sums checked
#> package 'jsonlite' successfully unpacked and MD5 sums checked
#> package 'curl' successfully unpacked and MD5 sums checked
#> package 'httr' successfully unpacked and MD5 sums checked
#> package 'ps' successfully unpacked and MD5 sums checked
#> package 'fastmap' successfully unpacked and MD5 sums checked
#> package 'digest' successfully unpacked and MD5 sums checked
#> package 'cachem' successfully unpacked and MD5 sums checked
#> package 'sass' successfully unpacked and MD5 sums checked
#> package 'tinytex' successfully unpacked and MD5 sums checked
#> package 'bslib' successfully unpacked and MD5 sums checked
#> package 'xfun' successfully unpacked and MD5 sums checked
#> package 'yaml' successfully unpacked and MD5 sums checked
#> package 'highr' successfully unpacked and MD5 sums checked
#> package 'evaluate' successfully unpacked and MD5 sums checked
#> package 'processx' successfully unpacked and MD5 sums checked
#> package 'rmarkdown' successfully unpacked and MD5 sums checked
#> package 'knitr' successfully unpacked and MD5 sums checked
#> package 'fs' successfully unpacked and MD5 sums checked
#> package 'callr' successfully unpacked and MD5 sums checked
#> package 'crayon' successfully unpacked and MD5 sums checked
#> package 'bit' successfully unpacked and MD5 sums checked
#> package 'vroom' successfully unpacked and MD5 sums checked
#> package 'tidyr' successfully unpacked and MD5 sums checked
#> package 'broom' successfully unpacked and MD5 sums checked
#> package 'timechange' successfully unpacked and MD5 sums checked
#> package 'readr' successfully unpacked and MD5 sums checked
#> package 'forcats' successfully unpacked and MD5 sums checked
#> package 'gargle' successfully unpacked and MD5 sums checked
#> package 'colorspace' successfully unpacked and MD5 sums checked
#> package 'isoband' successfully unpacked and MD5 sums checked
#> package 'data.table' successfully unpacked and MD5 sums checked
#> package 'pkgload' successfully unpacked and MD5 sums checked
#> package 'testthat' successfully unpacked and MD5 sums checked
#> package 'gtools' successfully unpacked and MD5 sums checked
#> package 'formatR' successfully unpacked and MD5 sums checked
#> package 'readxl' successfully unpacked and MD5 sums checked
#> package 'ragg' successfully unpacked and MD5 sums checked
#> package 'modelr' successfully unpacked and MD5 sums checked
#> package 'lubridate' successfully unpacked and MD5 sums checked
#> package 'ggplot2' successfully unpacked and MD5 sums checked
#> package 'dtplyr' successfully unpacked and MD5 sums checked
#> package 'dbplyr' successfully unpacked and MD5 sums checked
#> package 'pbapply' successfully unpacked and MD5 sums checked
#> package 'matrixStats' successfully unpacked and MD5 sums checked
#> package 'tidyverse' successfully unpacked and MD5 sums checked
#> package 'e1071' successfully unpacked and MD5 sums checked
#> 
#> The downloaded binary packages are in
#>  C:\Users\ypan\AppData\Local\Temp\RtmpiwZWXh\downloaded_packages
#>          checking for file 'C:\Users\ypan\AppData\Local\Temp\RtmpiwZWXh\remotes231852056eda\YuePan027-mcQTL-4bb8ffe/DESCRIPTION' ...  ✔  checking for file 'C:\Users\ypan\AppData\Local\Temp\RtmpiwZWXh\remotes231852056eda\YuePan027-mcQTL-4bb8ffe/DESCRIPTION'
#>       ─  preparing 'mcQTL':
#>    checking DESCRIPTION meta-information ...     checking DESCRIPTION meta-information ...   ✔  checking DESCRIPTION meta-information
#>       ─  checking for LF line-endings in source and make files and shell scripts
#> ─  checking for empty or unneeded directories
#>       ─  building 'mcQTL_0.1.0.tar.gz'
#>      
#> 
library(mcQTL)
library(ggplot2)
library(reshape2)
```

# Quick start

## Cell-type proportion deconvolution

This step is used to obtain cell-type proportion.

A `SummarizedExperiment` object with bulk protein or gene expression
contained in `counts` slot, and a “signature matrix” which serves as a
reference of known cellular signatures contained as an element in
`metadata` slot is required as input file. Note that the proteins or
genes in signature matrix can be different from that in assay, but only
common proteins or genes will be used to do deconvolution.

In this current version, only `CIBERSORT` and `nnls` are supported as
the deconvolution methods.

``` r
se <- SummarizedExperiment(assays = list(counts = mcQTL::protein_data),
                           rowData = mcQTL::anno_protein)
metadata(se) <- list(sig_matrix = mcQTL::ref_data)
se <- deconv(se, "cibersort")
```

This step might take a few minutes if there are many proteins or genes
in the signature matrix. The cell-type proportion estimates for each
sample will be stored as an element (`prop`) in `metadata` slot.

``` r
head(se@metadata$prop)
#>          CellType_1 CellType_2 CellType_3 CellType_4
#> Sample_1  0.2325557 0.06702399  0.5997269 0.10069340
#> Sample_2  0.1756690 0.04542362  0.6332230 0.14568439
#> Sample_3  0.1910644 0.08487123  0.6093489 0.11471552
#> Sample_4  0.2513749 0.07693023  0.5695605 0.10213434
#> Sample_5  0.2137109 0.10576155  0.5486246 0.13190293
#> Sample_6  0.3227685 0.05524151  0.5693060 0.05268406

ggplot(data.frame(reshape2::melt(se@metadata$prop)), 
       aes(x = Var2, y = value, fill = Var2)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1,
                                           dodge.width = 0.7),
           aes(fill = Var2, color = Var2),
           pch = 21, alpha = 0.5) +
  geom_boxplot(lwd=0.7, outlier.shape = NA) +
  theme_classic() +
  xlab("Cell type") + ylab("Estimated proportion")
```

![](README_files/figure-gfm/decov2-1.png)<!-- -->

Alternatively, if there are cell-type proportion estimates results
generated using other methods or obtained from other sources, just save
that as an element (`prop`) in `metadata` slot and this deconvolution
step can be skipped. Note that the samples in the cell-type proportion
estimates must match the samples from bulk protein/gene expression data.

## Feature filtering

The feature filtering can be applied at both proteins/genes and SNPs.
This step is optional but highly recommended to filter out some features
that are not very informative or do not make much sense biologically.
Note that this function is required to run even no filtering is expected
to be done (just set `filter_method = "null"`) to obtain a consistent
object format for downstream analysis.

To apply feature filtering, annotation files for protein/gene and SNPs
are required. The annotation file for proteins/genes should be stored in
`rowData()`, where each row corresponds to a protein/gene with it’s
symbol as row names. The first column should be a character vector
indicating which chromosome each protein or gene is on. In addition, it
should contain at least a “Start” column with numeric values indicating
the start position on that chromosome, a “End” column with numeric
values indicating the end position on that chromosome and a “Symbol”
column as a unique name for each protein or gene.

``` r
head(rowData(se))
#> DataFrame with 6 rows and 4 columns
#>                   Chr     Start       End      Symbol
#>           <character> <integer> <integer> <character>
#> Protein_1          15  67202823  67254631   Protein_1
#> Protein_2           6  44300549  44313323   Protein_2
#> Protein_3           7 122076491 122133726   Protein_3
#> Protein_4          16   8735739   8781427   Protein_4
#> Protein_5           9 104784317 104903679   Protein_5
#> Protein_6           9 137007931 137028140   Protein_6
```

The information from genetic variants should be stored in a P (the
number of SNP) by N (the number of samples, should match the sample in
`counts` slot) matrix contained as an element (`SNP_data`) in `metadata`
slot. Each matrix entry corresponds to the genotype group indicator (0
for 0/0, 1 for 0/1 and 2 for 1/1) for a sample at a genetic location.
The annotations of these SNP should be stored as an element (`anno_SNP`)
in `metadata` slot. It should include at least the following columns:
(1) “CHROM” (which chromosome the SNP is on); (2) “POS” (position of
that SNP) and (3) “ID” (a unique identifier for each SNP, usually a
combination of chromosome and its position).

The example SNP data provided here were restricted to chromosome 9 only.
In practice, the SNPs may from multiple or even all chromosomes.

``` r
se@metadata$SNP_data <- mcQTL::SNP_data
se@metadata$anno_SNP <- mcQTL::anno_SNP
head(se@metadata$anno_SNP)
#>        CHROM       POS          ID
#> 237392     9 104596634 9:104596634
#> 106390     9  28487163  9:28487163
#> 304108     9 126307371 9:126307371
#> 295846     9 122787821 9:122787821
#> 126055     9  33975396  9:33975396
#> 342900     9 140300675 9:140300675
```

For filtering at protein or gene level, only those symbols contained in
`target_protein` argument will be kept for csQTL analysis in the next
step. By default, all proteins or genes will be used.

For filtering at SNP level, there are three options: (1) only those
symbols contained in `target_SNP` argument will be kept and if not
provided, all SNPs will be used for further filtering; (2) filter out
the SNPs that have minor allele frequency below the threshold defined by
`filter_allele` argument (`filter_method = "allele"`) and (3) restrict
to cis-regulatory variants, i.e. the SNPs up to 1 Mb proximal to the
start of the gene (`filter_method = "distance"`).

The results after filtering will be stored as an element
(`choose_SNP_list`) in `metadata` slot. It is a list with the length of
the number of proteins for downstream analysis. Each element stores the
index of SNPs to be tested for corresponding protein. The proteins with
no SNPs correspond to it will be removed from the returned list.

To simplify the analysis, we only kept 10 targeted proteins from
chromosome 9 as an example.

``` r
target_protein <- rowData(se)[rowData(se)$Chr == 9,][1:10, "Symbol"]
se <- feature_filter(se, target_protein = target_protein, 
                     filter_method = c("allele", "distance"), 
                     filter_allele = 0.25,
                     filter_geno = 0.05,
                     ref_position = "TSS")           
#> Filter SNP based on distance for protein Protein_261
#> Filter SNP based on distance for protein Protein_241
#> Filter SNP based on distance for protein Protein_238
#> Filter SNP based on distance for protein Protein_131
#> Filter SNP based on distance for protein Protein_93
#> Filter SNP based on distance for protein Protein_88
#> Filter SNP based on distance for protein Protein_79
#> Filter SNP based on distance for protein Protein_6
#> Filter SNP based on distance for protein Protein_5
#> Filter SNP based on distance for protein Protein_283
```

In this example, only 20 SNPs are kept for the first target protein and
only 20 SNPs are kept for the second target protein.

``` r
unlist(lapply(se@metadata$choose_SNP_list, length))
#>   Protein_5   Protein_6  Protein_79  Protein_88  Protein_93 Protein_238 
#>          20          20           6           9           9           9 
#> Protein_241 Protein_261 Protein_283 
#>          23          12          15
```

## csQTL analysis

In this step, the `TOAST` method is implemented for cell-type-specific
differential expression analysis based on samples’ genotype.

The result will be stored as an element (`TOAST_output`) in `metadata`
slot. It is a list with the same length as tested proteins or genes
where each element consists of a table including protein or gene symbol,
SNP ID and p-values from each cell type. A significant p-value indicates
that the protein or gene expression is different among the sample from
different genotype groups.

``` r
system.time(se <- csQTL(se))
#> csQTL test for protein Protein_5 
#> csQTL test for protein Protein_6 
#> csQTL test for protein Protein_79 
#> csQTL test for protein Protein_88 
#> csQTL test for protein Protein_93 
#> csQTL test for protein Protein_238 
#> csQTL test for protein Protein_241 
#> csQTL test for protein Protein_261 
#> csQTL test for protein Protein_283
#>    user  system elapsed 
#>    3.55    1.28  457.83
```

We can check the results from csQTL analysis for the first target
protein by calling:

``` r
head(se@metadata$TOAST_output[[1]])
#>     protein         SNP CellType_1 CellType_2 CellType_3 CellType_4
#> 1 Protein_5 9:104702846  0.6687429  0.6076771  0.7661268 0.85811387
#> 2 Protein_5 9:104294530  0.9907045  0.8567053  0.8213907 0.47889613
#> 3 Protein_5 9:104868031  0.2276643  0.6469674  0.5206621 0.04218665
#> 4 Protein_5 9:105077464  0.7764349  0.9948742  0.8729841 0.27895985
#> 5 Protein_5 9:105655745  0.5509189  0.8004329  0.5812060 0.01227444
#> 6 Protein_5 9:105510462  0.3154370  0.2156796  0.3124176 0.01278241
```

## TCA deconvolution (optional)

The cell-type-specific gene expression per bulk sample can also be
estimated using `TCA` deconvolution method given cellular composition
(stored in `prop` in `metadata` slot). The output will be stored as an
element (`TCA_deconv`) in `metadata` slot. It is a list with the length
of the number of cell types (same as cell types in `prop` in `metadata`
slot). Each element stores a deconvoluted protein expression per bulk
sample. Below is an example to check the deconvoluted cellular
expression for the first cell type (restricted to first 5 proteins and
first 5 samples):

``` r
se <- TCA_deconv(se)
#> INFO [2023-02-27 15:22:29] Validating input...
#> INFO [2023-02-27 15:22:29] Starting tensor for estimating Z...
#> INFO [2023-02-27 15:22:29] Estimate tensor...
#> INFO [2023-02-27 15:22:32] Finished estimating tensor.
se@metadata$TCA_deconv[[1]][1:5,1:5]
#>           Sample_1 Sample_2 Sample_3 Sample_4 Sample_5
#> Protein_1 15.76586 16.08454 15.82306 15.76538 15.63249
#> Protein_2 16.34240 16.34232 16.34242 16.34250 16.34245
#> Protein_3 22.28076 21.47328 21.38018 22.64709 21.89433
#> Protein_4 19.48600 19.38766 19.54709 19.64752 19.61506
#> Protein_5 19.51987 19.36304 17.68337 19.18520 19.01969
```
