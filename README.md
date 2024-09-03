# AnnotateThat

An R package designed to create IonReporter custom annotation files for specified gene panels.

Supported annotation sources:
- gnomAD
- ClinVar
- InterVar 
- AlphaMissense
- MutationTaster
- Revel

## Installation from GitHub

Use with caution as the latest version may be unstable, although typical Bioconductor checks are executed before each push.

```r
if (!requireNamespace('devtools',quietly=TRUE))
    install.packages('devtools')

library(devtools)
install_github('BiobankCy/AnnotateThat')
```

## Basic use

```r
library(AnnotateThat)
gns <- c('SRY', 'AMELY') # a vector of genes
annotate(gns, annotators = 'revel') # create a custom REVEL vcf file 
annotate(gns, annotators = 'alphamissense') # create a custom AlphaMissense vcf file 
annotate(gns, annotators = 'intervar', databases = 'all', type = 'both') # create a custom InterVar prediction vcf file 
annotate(gns, annotators = 'gnomad_af', databases = 'gnomad_auto', type = 'genome') # create a custom gnomAD allele frequencies (AF) vcf file auto-downloaded from gnomAD genome files 
annotate(gns, annotators = 'mutation_taster', databases = 'all', type = 'both') # create a custom MutationTaster predictions vcf file 
annotate(gns, annotators = 'all', databases = 'all', type = 'both') # create all the above custom vcf files
```

Argument `databases` defines the online variant sources used to collect genetic variations for the gene vector provided. It is taken into consideration when `intervar`, `gnomad_af`, `mutation_taster` or `all` annotators are selected.

Argument `type` defines whether `exome`, `genome` or `both` genomics data types should be retrieved from gnomAD. It is taken into consideration in parallel with `databases` argument.

## Reference genome version

AnnotateThat annotation is by default built around **GRCh38**. For **GRCh37** variants, set `liftover = TRUE`.

## System requirements

Enough storage space should be available for downloading big gnomAD vcf files. For more details regarding file size please refer to respective [download page](https://gnomad.broadinstitute.org/downloads).