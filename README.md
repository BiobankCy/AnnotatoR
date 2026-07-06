# AnnotatoR

<!-- badges: start -->
  [![License](https://img.shields.io/badge/license-GPL--3.0-orange)](https://github.com/BiobankCy/AnnotatoR/blob/main/LICENSE)
  [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21098338.svg)](https://doi.org/10.5281/zenodo.21098338)
  ![GitHub repo size](https://img.shields.io/github/repo-size/BiobankCy/AnnotatoR)
  ![GitHub issues](https://img.shields.io/github/issues/BiobankCy/AnnotatoR)
  ![contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)
<!-- badges: end -->

An R package designed to create general purpose IonReporter-compatible annotation files starting from a set of gene names.

## Installation from GitHub

```r
if (!requireNamespace('pak',quietly=TRUE))
    install.packages('pak')

pak::pkg_install('github::BiobankCy/AnnotatoR')
```

## Basic use

`annotate()` is the main workhorse of the package. Through this wrapper, the user can provide a set of gene symbols to download annotation from:
* [REVEL](https://sites.google.com/site/revelgenomics/)
* [AlphaMissense](https://alphamissense.hegelab.org/)
* [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/)
* [gnomAD](https://gnomad.broadinstitute.org/)
* [InterVar](https://wintervar.wglab.org/)
* [MutationTaster](https://genecascade.org/)

for gene variants retrieved through:
* [gnomAD](https://gnomad.broadinstitute.org/)
* [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/)
* [LOVD3](https://www.lovd.nl/3.0/home)

### Arguments

* `gns` : a set of gene names
* `annotators` : defines the annotator source(s) to use (*revel*, *alphamissense*, *clinvar_sig*, *gnomad_af*, *intervar*, 
    *mutation_taster*, *all*)
* `databases` : defines the online, publicly available databases to retrieve variants from (*gnomad_man*, *gnomad_auto*, *clinvar*, *lovd3*, *all*). The argument is taken into consideration when *intervar*, *gnomad_af*, *mutation_taster* or *all* annotators are selected. *gnomad_man* expects manually downloaded gnomAD csv files placed in *gnomad/manual* subdirectory of the directory set in *path* argument (see below).
* `panelName` : defines the path and prefix of the output file
* `path` : defines the path were database files are stored (default: *./dbs*)
* `type` : defines which genomic data types to retrieve from gnomAD (*exome*, *genome* or *both*). It is taken into consideration in parallel with *databases* argument.
* `liftover` : defines whether variants should be lifted over from GRCh38 to GRCh37 (default: FALSE)
* `saveRaw` : defines whether retrieved variants should be saved without any annotation. Applicable only when *intervar* or *all* annotator is selected.

### Examples

The basic input is a set of genes
```r
gns <- c('SRY', 'AMELY')
```

Retrieve REVEL predictions for all variants in `gns` as listed in AlphaMissense pre-computed file
```r
annotate(gns, annotators = 'revel')
```

Retrieve AlphaMissense predictions for all variants in `gns` as listed in AlphaMissense pre-computed file
```r
annotate(gns, annotators = 'alphamissense') # create a custom AlphaMissense vcf file 
```

Retrieve InterVar-calculated ACMG predictions for all variants in `gns` as retrieved from `all` databases
```r
annotate(gns, annotators = 'intervar', databases = 'all', type = 'both') # create a custom InterVar prediction vcf file 
```

Retrieve InterVar-calculated ACMG predictions for all variants in `gns` as retrieved from `all` databases
```r
annotate(gns, annotators = 'intervar', databases = 'all', type = 'both') # create a custom InterVar prediction vcf file 
```

Retrieve gnomAD allele frequencies for all variants in `gns` as retrieved from `gnomad` database
```r
annotate(gns, annotators = 'gnomad_af', databases = 'gnomad_auto', type = 'genome')
```

Retrieve gnomAD allele frequencies for all variants in `gns` as retrieved from manually fetched gnomAD csv files stored in `dbs/gnomad/manual`
```r
annotate(gns, annotators = 'gnomad_af', databases = 'gnomad_man')
```

Retrieve MutationTaster2021 (GRCh37) or MutationTaster2025 (GRCh38) allele frequencies for all variants in `gns` as retrieved from `all` databases. `both` exome and genome based gnomAD files are used.
```r
annotate(gns, annotators = 'mutation_taster', databases = 'all', type = 'both')
```

Retrieve annotation from all supported annotators and databases. 
```r
annotate(gns, annotators = 'all', databases = 'all', type = 'both')
```

## Required storage space

AnnotatoR can use **extensive disk space**, depending on its use. For this reason, user is asked **each time** a file is about to be downloaded.

### gnomad_man vs gnomad_auto

Whenever gnomAD-related annotation is required, e.g. to retrieve allele frequencies or ClinVar significance, data can be retrieved from either gnomAD download website (`gnomad_auto` databases option), or from manually downloaded csv files placed in `gnomad/manual` subdirectory of the `path` directory (`gnomad_man` databases option). The former option requires significant empty space on the disk and should be avoided, when possible. For small gene sets, the latter option is more viable and suggested. By default, the package searches for manually downloaded files and if not present asks to download gnomAD chromosome files instead.

For more details regarding gnomAD file size please refer to the respective [download page](https://gnomad.broadinstitute.org/downloads).

## Reference genome version

AnnotatoR annotation is by default built around **GRCh38**. For **GRCh37** variants, set `liftover = TRUE`.

## Funding

This project has received funding from the European Union's Horizon 2020 research and innovation programme under grant agreement No 857122 ([CY-BIOBANK, Biobanking and the Cyprus Human Genome Project](https://cordis.europa.eu/project/id/857122)).