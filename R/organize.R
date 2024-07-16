`%!in%` <- Negate(`%in%`)

#' Contruct AlphaMissense annotation file.
#'
#' Download and organize AlphaMissense predictions for specific genes.
#'
#' @param gns Gene names character vector.
#' @param path Where data were downloaded
#' @return vcf data.frame
#' @export
constructAM <- function(gns, path){
	path <- file.path(path, 'alphamissense')
	suppressWarnings(dir.create(path, recursive = TRUE))
	message('Downloading AlphaMissense')
	download_am(path)

	edb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
	granges <- ensembldb::genes(edb)
	tranges <- ensembldb::transcripts(edb)
	gns_id <- unique(as.data.frame(granges[GenomicRanges::elementMetadata(granges)$symbol %in% gns])[,'gene_id'])
	gns_id <- gns_id[grep('ENSG', gns_id)]
	map <- as.data.frame(tranges[GenomicRanges::elementMetadata(tranges)$gene_id %in% gns_id])
	map$seqnames <- as.character(map$seqnames)

	message('Subsetting AlphaMissense')
	system(paste0("gzip -cdk ", file.path(path, "AlphaMissense_hg38.tsv.gz"),
		" | LC_ALL=C grep -i -E '", paste0(map$tx_id, collapse = '|'), 
		"' > ", file.path(path, "alphaMissense_sel.txt")))
	system(paste0("gzip -cdk ", file.path(path, "AlphaMissense_hg38.tsv.gz"),
		" | head -n 4 | tail -n 1 > ", file.path(path, "alphaMissense_header.txt")))

	am <- read.delim(file.path(path, 'alphaMissense_sel.txt'), header = FALSE)
	colnames(am) <- read.delim(file.path(path, 'alphaMissense_header.txt'), header = FALSE)
	am$transcript_id <- gsub('\\..*', '', am$transcript_id)
	am <- merge(am, map, by.x = 'transcript_id', by.y = 'tx_id')
	# Round-up am pathogenicity score (for easier filter chain creation)
	am$am_pathogenicity <- round(am$am_pathogenicity,2)

	am_vcf <- data.frame(CHR = am[,'#CHROM'], POS = am$POS,
		ID = paste0(stringr::str_to_title(am$am_class), '(', am$am_pathogenicity, ')'),
		REF = am$REF, ALT = am$ALT, QUAL = '.', FILTER = 'PASS', INFO = '.',
		FORMAT = '.', Sample = '.'
	)
	return(am_vcf)
}

#' Contruct REVEL annotation file.
#'
#' Download and organize REVEL meta-predictions for specific genes.
#'
#' @param gns Gene names character vector.
#' @param path Where data were downloaded
#' @return vcf data.frame
#' @export
constructREVEL <- function(gns, path){
	path <- file.path(path, 'revel')
	suppressWarnings(dir.create(path, recursive = TRUE))

	# Gene coordinates (hg19 for revel 1.3 file names)
	edb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
	granges <- ensembldb::genes(edb)
	map <- as.data.frame(
		granges[GenomicRanges::elementMetadata(granges)$symbol %in% gns])
	map <- map[grep('ENSG', map$gene_id), ]
	map$seqnames <- as.character(map$seqnames)
	map <- split(map, f = map$seqnames)

	# Download revel data files for respective chromosomes
	message('Fetching and Subsetting REVEL files.')
	revel_data <- list()
	for(i in names(map)){
		if(nchar(i) == 1 & i %!in% c('X', 'Y')){
			ii <- paste0('0', i)
		} else if(nchar(i) == 1 & i %in% c('X', 'Y')){
			ii <- paste0('_', i)
		} else {
			ii <- i
		}
		name <- download_revel(ii, path)
		
		# Find the correct gene segment(s)
		segs <- do.call('rbind', lapply(strsplit(
			list.files(name), split = '_'), function(x){
				as.data.frame(t(as.numeric(gsub('.csv', '', x[(length(x)-1):length(x)]))))
			}
		))
		colnames(segs) <- c('start', 'end')
		rownames(segs) <- NULL

		# Gather data per gene
		for(k in rownames(map[[i]])){
			revel <- read.csv(
				list.files(name, full.names = TRUE)[which(
					segs$start < map[[i]][k, 'start'] & 
					segs$end > map[[i]][k, 'end'])]
			)
			revel_vcf <- data.frame(CHR = paste0('chr', revel$chr),
				POS = revel$grch38_pos, ID = round(revel$REVEL, 2), 
				REF = revel$ref, ALT = revel$alt, QUAL = '.',
				FILTER = 'PASS', INFO = '.', FORMAT = '.', Sample = '.')
			# # Maintain variants within panel genes start/end +/- 1000 nts
			# revel_data[[k]] <- revel_vcf[which(revel_vcf$POS > map[[i]][k, 'start'] - 1000 & 
			# 	revel_vcf$POS > map[[i]][k, 'end'] + 1000),]
			revel_data[[i]] <- revel_vcf

		}
	}
	revel_data <- do.call('rbind', revel_data)
	rownames(revel_data) <- NULL
	revel_data <- unique(revel_data)
	return(revel_data)
}


#' Contruct ClinVar annotation file.
#'
#' Download and organize ClinVar significance classification with
#' review status level for specific genes.
#'
#' @param gns Gene names character vector.
#' @param path Where data were downloaded
#' @return vcf data.frame
#' @export
constructClinVar <- function(gns, path){
	path <- file.path(path, 'clinvar')
	suppressWarnings(dir.create(path, recursive = TRUE))

	message('Fetching clinvar variants')
	download_clinvar(path)
	message('Subsetting clinvar variants')
	file.gz <- file.path(path, 'clinvar.vcf.gz')
	file.gz.tbi <- paste(file.gz, ".tbi", sep="")
	if(!(file.exists(file.gz.tbi)))
	    Rsamtools::indexTabix(file.gz, format="vcf")
	# CLNREVSTAT: review status for the aggregate germline classification
	params <- VariantAnnotation::ScanVcfParam(info = c('GENEINFO', 'CLNSIG', 'CLNREVSTAT'))
	vcf <- VariantAnnotation::readVcf(Rsamtools::TabixFile(file.gz), 'hg38', params)
	rg <- SummarizedExperiment::rowRanges(vcf)
	inf <- VariantAnnotation::info(vcf)
	clinvar_vcf <- data.frame(
		CHR = as.character(GenomicRanges::seqnames(rg)), 
		POS = as.data.frame(rg@ranges)$start,
		REF = as.character(S4Vectors::DataFrame(rg)$REF), 
		ALT = as.character(unlist(S4Vectors::DataFrame(rg)$ALT)),
		GENEINFO = as.character(inf$GENEINFO),
		SIG = as.character(unlist(lapply(
			S4Vectors::DataFrame(inf)$CLNSIG, paste, collapse = ''))),
		REVSTAT= as.character(unlist(lapply(
			S4Vectors::DataFrame(inf)$CLNREVSTAT, paste, collapse = ',')))
	)
	# Genes of interest & at least one star
	clinvar_vcf <- merge(
		clinvar_vcf[grep(paste(gns, collapse = '|'), clinvar_vcf$GENEINFO), ],
		revstatus_map[which(revstatus_map$stars != 'none'), ], 
		by = 'REVSTAT'
	)
	clinvar_vcf <- data.frame(CHR = paste0('chr', clinvar_vcf$CHR),
		POS = clinvar_vcf$POS, ID = paste0(clinvar_vcf$SIG, '(', clinvar_vcf$stars, ')'), 
		REF = clinvar_vcf$REF, ALT = clinvar_vcf$ALT, QUAL = '.',
		FILTER = 'PASS', INFO = '.', FORMAT = '.', Sample = '.')
	return(clinvar_vcf)
}

#' Collect known variants.
#'
#' Collect known variants for specific genes from public resources
#' (gnomAD, ClinVar, LOVD3).
#' 
#' @param gns Gene names character vector.
#' @param databases Databases to retrieve variants from. 
#' @param path Where data were downloaded
#' @param type Type of gnomAD data to download.
#' @return vcf data.frame
#' @export
collectVars <- function(gns, databases = c('gnomad', 'clinvar', 'lovd3', 'all'), 
	path, type = c('exomes', 'genomes', 'both')){

	# Get gene data
	message('Preparing data')
	edb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
	granges <- ensembldb::genes(edb)
	map <- as.data.frame(granges[GenomicRanges::elementMetadata(granges)$symbol %in% gns])
	map <- map[grep('ENSG', map$gene_id), ]
	map$seqnames <- as.character(map$seqnames)
	map <- split(map, f = map$seqnames)

	# gnomad
	pathTmp <- file.path(path, 'gnomad')
	suppressWarnings(dir.create(pathTmp, recursive = TRUE))
	message('Fetching gnomad variants')
	gnomad_vcf <- do.call('rbind', lapply(map, function(x){
		switch(type,
			exomes = {
				download_gnomad(paste0('chr', unique(x$seqnames)), type, pathTmp)
				out <- list()
				for(k in type) {
					message('Subsetting gnomAD ', k)
					file.gz <- file.path(pathTmp, paste0('gnomad.', k, '.v4.0.sites.chr', 
						unique(x$seqnames), '.vcf.bgz'))
					file.gz.tbi <- paste(file.gz, ".tbi", sep="")
					gr <- GenomicRanges::GRanges(paste0('chr', x$seqnames), IRanges::IRanges(x$start, x$end))
					params <- VariantAnnotation::ScanVcfParam(which = gr)
					vcf <- VariantAnnotation::readVcf(Rsamtools::TabixFile(file.gz), 'hg38', params)
					rg <- SummarizedExperiment::rowRanges(vcf)
					out[[k]] <- data.frame(
						CHR = as.character(GenomicRanges::seqnames(rg)), 
						POS = as.data.frame(rg@ranges)$start, 
						REF = as.character(S4Vectors::DataFrame(rg)$REF), 
						ALT = as.character(unlist(S4Vectors::DataFrame(rg)$ALT))
					)
				}
				return(do.call('rbind', out))
			},
			genomes = {
				download_gnomad(paste0('chr', unique(x$seqnames)), type, pathTmp)
				out <- list()
				for(k in type) {
					message('Subsetting gnomAD ', k)
					file.gz <- file.path(pathTmp, paste0('gnomad.', k, '.v4.0.sites.chr', 
						unique(x$seqnames), '.vcf.bgz'))
					file.gz.tbi <- paste(file.gz, ".tbi", sep="")
					gr <- GenomicRanges::GRanges(paste0('chr', x$seqnames), IRanges::IRanges(x$start, x$end))
					params <- VariantAnnotation::ScanVcfParam(which = gr)
					vcf <- VariantAnnotation::readVcf(Rsamtools::TabixFile(file.gz), 'hg38', params)
					rg <- SummarizedExperiment::rowRanges(vcf)
					out[[k]] <- data.frame(
						CHR = as.character(GenomicRanges::seqnames(rg)), 
						POS = as.data.frame(rg@ranges)$start, 
						REF = as.character(S4Vectors::DataFrame(rg)$REF), 
						ALT = as.character(unlist(S4Vectors::DataFrame(rg)$ALT))
					)
				}
				return(do.call('rbind', out))
			},
			both = {
				types <- c('exomes', 'genomes')
				download_gnomad(paste0('chr', unique(x$seqnames)), types, pathTmp)
				out <- list()
				for(k in types) {
					message('Subsetting gnomAD ', k)
					file.gz <- file.path(pathTmp, paste0('gnomad.', k, '.v4.0.sites.chr', 
						unique(x$seqnames), '.vcf.bgz'))
					file.gz.tbi <- paste(file.gz, ".tbi", sep="")
					gr <- GenomicRanges::GRanges(paste0('chr', x$seqnames), IRanges::IRanges(x$start, x$end))
					params <- VariantAnnotation::ScanVcfParam(which = gr)
					vcf <- VariantAnnotation::readVcf(Rsamtools::TabixFile(file.gz), 'hg38', params)
					rg <- SummarizedExperiment::rowRanges(vcf)
					out[[k]] <- data.frame(
						CHR = as.character(GenomicRanges::seqnames(rg)), 
						POS = as.data.frame(rg@ranges)$start, 
						REF = as.character(S4Vectors::DataFrame(rg)$REF), 
						ALT = as.character(unlist(S4Vectors::DataFrame(rg)$ALT))
					)
				}
				return(do.call('rbind', out))
			}

		)
	}))
	rownames(gnomad_vcf) <- NULL

	# clinvar
	pathTmp <- file.path(path, 'clinvar')
	suppressWarnings(dir.create(pathTmp, recursive = TRUE))
	message('Fetching clinvar variants')
	download_clinvar(pathTmp)
	message('Subsetting clinvar variants')
	file.gz <- file.path(pathTmp, 'clinvar.vcf.gz')
	file.gz.tbi <- paste(file.gz, ".tbi", sep="")
	if(!(file.exists(file.gz.tbi)))
	    Rsamtools::indexTabix(file.gz, format="vcf")
	params <- VariantAnnotation::ScanVcfParam(info = 'GENEINFO')
	vcf <- VariantAnnotation::readVcf(Rsamtools::TabixFile(file.gz), 'hg38', params)
	rg <- SummarizedExperiment::rowRanges(vcf)
	inf <- VariantAnnotation::info(vcf)
	clinvar_vcf <- data.frame(
		CHR = as.character(GenomicRanges::seqnames(rg)), 
		POS = as.data.frame(rg@ranges)$start, 
		REF = as.character(S4Vectors::DataFrame(rg)$REF), 
		ALT = as.character(unlist(S4Vectors::DataFrame(rg)$ALT)),
		GENEINFO = as.character(inf$GENEINFO)
	)
	# system(paste0('bcftools query -f "%CHROM %POS %REF %ALT %GENEINFO\n" ', file.path(pathTmp, 'clinvar.vcf.gz'), 
	# 	' > ', file.path(pathTmp, 'clinvar_tmp.vcf')))
	# clinvar_vcf <- read.table(file.path(pathTmp, 'clinvar_tmp.vcf'))
	# colnames(clinvar_vcf) <- c('CHR', 'POS', 'REF', 'ALT', 'GENEINFO')
	clinvar_vcf <- clinvar_vcf[grep(paste(gns, collapse = '|'), clinvar_vcf$GENEINFO), 1:4]

	# lovd3
	pathTmp <- file.path(path, 'lovd3')
	suppressWarnings(dir.create(pathTmp, recursive = TRUE))
	message('Retrieving LOVD3 variants')
	download_lovd3(gns, pathTmp)
	lovd3_vcf <- list()
	for(i in list.files(pathTmp, full.names = TRUE)){
		tmp <- readLines(i)
		gene <- gsub('LOVD_full_download_|_.*-.*', '', basename(i))

		# Gene without public variants?
		if(length(grep('^Error', tmp))){
			message('No public variants for ', gene)
			lovd3_vcf[[gene]] <- NA
		} else {
			message('Now processing ', gene)
			# Genomic data
			start <- grep('\\#\\# Variants_On_Genome', tmp)
			stop <- grep('\\#\\# Variants_On_Transcripts', tmp) -4
			df_genom <- read.delim(i, nrows = stop - (start + 2), skip = start + 3, header = FALSE)
			colnames(df_genom) <- gsub('\\"\\{\\{|\\}\\}\\"', '', strsplit(tmp[start + 3], split = '\t')[[1]])

			# Transcript data
			start <- grep('\\#\\# Variants_On_Transcripts', tmp)
			stop <- grep('\\#\\# Screenings_To_Variants', tmp) -4
			df_trans <- read.delim(i, nrows = stop - (start + 3), skip = start + 4, header = FALSE)
			colnames(df_trans) <- gsub('\\"\\{\\{|\\}\\}\\"', '', strsplit(tmp[start + 4], split = '\t')[[1]])

			# Required columns
			lovd3_vcf[[gene]] <- unique(merge(df_genom, df_trans, by = 'id', 
				all.x = TRUE, all.y = TRUE)[,
				c('id', 'chromosome', 'VariantOnGenome/DNA/hg38', 
				'VariantOnGenome/Genetic_origin', 
				'VariantOnGenome/ClinicalClassification',
				'VariantOnGenome/Published_as', 'VariantOnTranscript/DNA')])
			lovd3_vcf[[gene]]$chromosome <- paste0('chr', lovd3_vcf[[gene]]$chromosome)
			lovd3_vcf[[gene]] <- cbind(gene = gene, lovd3_vcf[[gene]])
			lovd3_vcf[[gene]]$id <- NULL
			colnames(lovd3_vcf[[gene]]) <- gsub('/', '.', colnames(lovd3_vcf[[gene]]))
		}
	}
	lovd3_vcf[is.na(lovd3_vcf)] <- NULL
	if(length(lovd3_vcf)){
		lovd3_vcf <- do.call('rbind', lovd3_vcf)
		rownames(lovd3_vcf) <- NULL
		} else {
			lovd3_vcf <- data.frame(CHR = character(), POS = character(),
				REF = character(), ALT = character())
		}

	# Clean
	unlink(c(
		file.path(path, 'gnomad', 'tmp_gen.vcf'), 
		file.path(path, 'gnomad', 'tmp_ex.vcf'), 
		file.path(path, 'clinvar', 'clinvar_tmp.vcf'), 
		file.path(path, 'lovd3', '*')
	))

	out <- unique(rbind(gnomad_vcf, clinvar_vcf, lovd3_vcf))
	return(out)

}

