# AnnotateThat

An R package designed to create IonReporter custom annotation files for specified gene panels.

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
annotate(gns, annotators = 'intervar') # create a custom InterVar prediction vcf file 
annotate(gns, annotators = 'all') # create all the above custom vcf files
```

## System requirements

Enough storage space should be available for downloading big gnomAD vcf files. For more details regarding file size please refer to respective [download page](https://gnomad.broadinstitute.org/downloads).