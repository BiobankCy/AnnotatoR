`%!in%` <- Negate(`%in%`)
`%>%` <- dplyr::`%>%`

#' Contruct AlphaMissense annotation file.
#'
#' Download and organize AlphaMissense predictions for specific genes.
#'
#' @param gns Gene names character vector.
#' @param path Where data were downloaded
#' @param liftover Should variants be lifted over?
#' @return vcf data.frame
#' @export
constructAM <- function(gns, path, liftover){
	path <- file.path(path, 'alphamissense')
	suppressWarnings(dir.create(path, recursive = TRUE))
	message('Downloading AlphaMissense')
	download_am(path)

	message('Subsetting AlphaMissense')
	map <- map_fetch('EnsDb.Hsapiens.v86', gns, trans = TRUE)
	system(paste0("gzip -cdk ", file.path(path, "AlphaMissense_hg38.tsv.gz"),
		" | LC_ALL=C grep -i -E '", paste0(map$tx_id, collapse = '|'), 
		"' > ", file.path(path, "alphaMissense_sel.txt")))
	system(paste0("gzip -cdk ", file.path(path, "AlphaMissense_hg38.tsv.gz"),
		" | head -n 4 | tail -n 1 > ", file.path(path, "alphaMissense_header.txt")))

	am <- utils::read.delim(file.path(path, 'alphaMissense_sel.txt'), header = FALSE)
	colnames(am) <- utils::read.delim(file.path(path, 'alphaMissense_header.txt'), header = FALSE)
	am$transcript_id <- gsub('\\..*', '', am$transcript_id)
	am <- merge(am, map, by.x = 'transcript_id', by.y = 'tx_id')
	# Round-up am pathogenicity score (for easier filter chain creation)
	am$am_pathogenicity <- round(am$am_pathogenicity,2)

	am_vcf <- data.frame(CHR = am[,'#CHROM'], POS = am$POS,
		ID = paste0(stringr::str_to_title(am$am_class), '(', am$am_pathogenicity, ')'),
		REF = am$REF, ALT = am$ALT, QUAL = '.', FILTER = 'PASS', INFO = '.',
		FORMAT = '.', Sample = '.'
	)
	if(isTRUE(liftover)) am_vcf$POS <- as.data.frame(lift(am_vcf))$start
	return(am_vcf)
}

#' Contruct REVEL annotation file.
#'
#' Download and organize REVEL meta-predictions for specific genes.
#'
#' @param gns Gene names character vector.
#' @param path Where data were downloaded
#' @param liftover Should variants be lifted over?
#' @return vcf data.frame
#' @export
constructREVEL <- function(gns, path, liftover){
	path <- file.path(path, 'revel')
	suppressWarnings(dir.create(path, recursive = TRUE))

	# Gene coordinates (hg19 for revel 1.3 file names; see publication)
	map <- map_fetch('EnsDb.Hsapiens.v75', gns, trans = FALSE)
	map <- split(map, f = map$seqnames)

	# Download revel data files for respective chromosomes
	message('Fetching and Subsetting REVEL files')
	revel_data <- list()
	for(i in gsub('chr', '', names(map))){
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
		iii <- paste0('chr', i)
		for(k in rownames(map[[iii]])){
			# Control for a gene coordinates spanning over >1 revel files 
			start <- which(segs$start < map[[iii]][k, 'start'])
			end <- which(segs$end > map[[iii]][k, 'end'])
			revel_fls <- intersect(start, end)
			if(length(revel_fls)) {
				revel <- utils::read.csv(
					list.files(name, full.names = TRUE)[revel_fls]
				)
				} else {
					revel <- rbind( 
						utils::read.csv(list.files(name, full.names = TRUE)[tail(start, 1)]),
						utils::read.csv(list.files(name, full.names = TRUE)[head(end, 1)])
					)
				}
			# revel <- utils::read.csv(
			# 	list.files(name, full.names = TRUE)[which(
			# 		segs$start < map[[iii]][k, 'start'] & 
			# 		segs$end > map[[iii]][k, 'end'])]
			# )
			# Instead of applying a liftover, choose the hg19 revel column
			tmp <- ifelse(isTRUE(liftover), 'hg19_pos', 'grch38_pos')
			revel_vcf <- data.frame(CHR = paste0('chr', revel$chr),
				POS = revel[,tmp], ID = round(revel$REVEL, 2), 
				REF = revel$ref, ALT = revel$alt, QUAL = '.',
				FILTER = 'PASS', INFO = '.', FORMAT = '.', Sample = '.')
			# Maintain variants within panel genes
			revel_data[[iii]] <- revel_vcf[which(
				revel_vcf$POS > map[[iii]][k, 'start'] & 
				revel_vcf$POS < map[[iii]][k, 'end']),]
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
#' @param liftover Should variants be lifted over?
#' @return vcf data.frame
#' @export
constructClinVar_deprecated <- function(gns, path, liftover){
	path <- file.path(path, 'clinvar')
	suppressWarnings(dir.create(path, recursive = TRUE))
	message('Fetching ClinVar variants')
	download_clinvar(path)
	message('Subsetting ClinVar variants')
	file.gz <- file.path(path, 'clinvar.vcf.gz')
	file.gz.tbi <- paste(file.gz, ".tbi", sep="")
	if(!(file.exists(file.gz.tbi)))
	    Rsamtools::indexTabix(file.gz, format="vcf")
	# CLNREVSTAT: review status for the aggregate germline classification
	params <- VariantAnnotation::ScanVcfParam(info = c('GENEINFO', 'CLNSIG', 'CLNREVSTAT'))
	vcf <- VariantAnnotation::readVcf(Rsamtools::TabixFile(file.gz), 'hg38', params)
	# Genes of interest
	vcf_sel <- vcf[grep(paste0('^', gns, collapse = '|'), VariantAnnotation::info(vcf)$GENEINFO),]
	rg <- SummarizedExperiment::rowRanges(vcf_sel)
	inf <- VariantAnnotation::info(vcf_sel)
	clinvar_vcf <- data.frame(
		CHR = paste0('chr', as.character(GenomicRanges::seqnames(rg))), 
		POS = as.data.frame(rg@ranges)$start,
		REF = as.character(S4Vectors::DataFrame(rg)$REF), 
		ALT = as.character(unlist(S4Vectors::DataFrame(rg)$ALT)),
		GENEINFO = as.character(inf$GENEINFO),
		SIG = as.character(unlist(lapply(
			S4Vectors::DataFrame(inf)$CLNSIG, paste, collapse = ''))),
		REVSTAT= as.character(unlist(lapply(
			S4Vectors::DataFrame(inf)$CLNREVSTAT, paste, collapse = ',')))
	)
	# At least one star
	clinvar_vcf <- merge(
		clinvar_vcf[grep(paste0('^', gns, collapse = '|'), clinvar_vcf$GENEINFO), ],
		revstatus_map[which(revstatus_map$stars != 'none'), ], 
		by = 'REVSTAT'
	)
	# # Remove regions overlapping LOCs
	# clinvar_vcf <- clinvar_vcf[grep('LOC|DT', clinvar_vcf$GENEINFO, invert = TRUE), ]

	clinvar_vcf <- clinvar_vcf %>%
		dplyr::mutate(
		  	SIG = dplyr::recode(SIG, 
				'Conflicting_classifications_of_pathogenicity' = 'Conflicting',
				'Benign/Likely_benign' = 'BLB',
				'Likely_benign' = 'LB',
				'Benign' = 'B',
				'Uncertain_significance' = 'VUS',
				'Likely_pathogenic' = 'LP',
				'Pathogenic/Likely_pathogenic' = 'PLP',
				'Pathogenic' = 'P'
		  	)
		)
	clinvar_vcf <- data.frame(CHR = clinvar_vcf$CHR,
		POS = clinvar_vcf$POS, ID = paste0(clinvar_vcf$SIG, '(', clinvar_vcf$stars, ')'), 
		REF = clinvar_vcf$REF, ALT = clinvar_vcf$ALT, QUAL = '.',
		FILTER = 'PASS', INFO = '.', FORMAT = '.', Sample = '.')
	# # Remove off gene variants
	# map <- map_fetch('EnsDb.Hsapiens.v86', gns)
	# clinvar_vcf <- clinvar_vcf[which(clinvar_vcf$POS >= map$start & clinvar_vcf$POS <= map$end), ]
	if(isTRUE(liftover)) clinvar_vcf$POS <- as.data.frame(lift(clinvar_vcf))$start
	return(clinvar_vcf)
}

#' Retrieve ClinVar significance.
#'
#' A wrapper function to retrieve ClinVar variant significance from gnomAD.
#' Manually downloaded gnomAD files are required.
#'
#' @param gns Gene names character vector.
#' @param path Where data were downloaded
#' @param liftover Should variants be lifted over?
#' @return vcf data.frame
#' @export
constructClinVar <- function(gns, path, liftover){

	# Get gene data
	message('Preparing data')
	map <- map_fetch('EnsDb.Hsapiens.v86', gns, trans = FALSE)
	map <- split(map, f = map$seqnames)

	pathTmp <- file.path(path, 'gnomad/manual')
	if(!dir.exists(pathTmp)){
		stop('gnomAD files are missing. Manually downloaded gnomAD files are required.')
	} else if (dir.exists(pathTmp)) {
		gns_tmp <- do.call('rbind', map)$gene_id
		files <- list.files(pathTmp, pattern = paste('gnomAD_v4', gns_tmp, sep = '|'), 
			full.names = TRUE)
		if(length(grep(paste(gns_tmp, collapse = '|'), files)) == length(gns_tmp)){
			out <- gnomad_fetch_manual(map, path, liftover, sig = TRUE)
		} else {
			stop('gnomAD files are missing. Manually downloaded gnomAD files are required.')
		}
	}
	
	af_vcf <- data.frame(out, QUAL = '.', FILTER = 'PASS', INFO = '.', FORMAT = '.', Sample = '.')
	return(af_vcf)
}

#' Retrieve gnomAD allele frequencies.
#'
#' A wrapper function to retrieve population allele frequencies as recorded in gnomAD.
#' 
#' @param gns Gene names character vector.
#' @param path Where data were downloaded
#' @param type Type of gnomAD data to download.
#' @param liftover Should variants be lifted over?
#' @return vcf data.frame
#' @export
constructAF <- function(gns, path, databases = c('gnomad_man', 'gnomad_auto'), type = c('exomes', 'genomes', 'both'), liftover){

	# Get gene data
	message('Preparing data')
	map <- map_fetch('EnsDb.Hsapiens.v86', gns, trans = FALSE)
	map <- split(map, f = map$seqnames)

	if(databases %in% c('clinvar', 'lovd3', 'all')) databases <- 'gnomad_man'

	databases <- match.arg(databases)
	type <- match.arg(type)
	switch(databases,
		gnomad_man = {
			pathTmp <- file.path(path, 'gnomad/manual')
			if(!dir.exists(pathTmp)){
				message('gnomAD files are missing. Proceed to download')
				out <- gnomad_fetch_auto(map, type, path, liftover, af = TRUE)
			} else if (dir.exists(pathTmp)) {
				gns_tmp <- do.call('rbind', map)$gene_id
				files <- list.files(pathTmp, pattern = paste('gnomAD_v4', gns_tmp, sep = '|'), 
					full.names = TRUE)
				if(length(grep(paste(gns_tmp, collapse = '|'), files)) == length(gns_tmp)){
					out <- gnomad_fetch_manual(map, path, liftover, af = TRUE)
				} else {
					message('gnomAD files are missing. Proceed to download')
					out <- gnomad_fetch_auto(map, type, path, liftover, af = TRUE)
				}
			}
		},
		gnomad_auto = {
			out <- gnomad_fetch_auto(map, type, path, liftover, af = TRUE)
		}
	)
	af_vcf <- data.frame(out, QUAL = '.', FILTER = 'PASS', INFO = '.', FORMAT = '.', Sample = '.')
	return(af_vcf)
}

#' Collect known variants.
#'
#' A wrapper function to collect known variants for specific genes
#' from public resources (gnomAD, ClinVar, LOVD3).
#' 
#' @param gns Gene names character vector.
#' @param databases Databases to retrieve variants from. 
#' @param path Where data were downloaded
#' @param type Type of gnomAD data to download.
#' @param liftover Should variants be lifted over?
#' @return vcf data.frame
#' @export
collectVars <- function(gns, databases = c('gnomad_man', 'gnomad_auto', 'clinvar', 'lovd3', 'all'), 
	path, type = c('exomes', 'genomes', 'both'), liftover){

	# Get gene data
	message('Preparing data')
	map <- map_fetch('EnsDb.Hsapiens.v86', gns, trans = FALSE)
	map <- split(map, f = map$seqnames)

	databases <- match.arg(databases)
	type <- match.arg(type)
	switch(databases, 
		gnomad_man = {
			out <- gnomad_fetch_manual(map, path, liftover)
		},
		gnomad_auto = {
			out <- gnomad_fetch_auto(map, type, path, liftover)
		},
		clinvar = {
			progressr::with_progress(out <- clinvar_fetch(gns, path, liftover))
		},
		lovd3 = {
			progressr::with_progress(out <- lovd3_fetch(gns, path, liftover))
		},
		all = {
			pathTmp <- file.path(path, 'gnomad/manual')
			if(!dir.exists(pathTmp)){
				message('gnomAD files are missing. Proceed to download')
				a <- gnomad_fetch_auto(map, type, path, liftover)
			} else if (dir.exists(pathTmp)) {
				gns_tmp <- do.call('rbind', map)$gene_id
				files <- list.files(pathTmp, pattern = paste('gnomAD_v4', gns_tmp, sep = '|'), 
					full.names = TRUE)
				if(length(grep(paste(gns_tmp, collapse = '|'), files)) == length(gns_tmp)){
					a <- gnomad_fetch_manual(map, path, liftover)
				} else {
					message('gnomAD files are missing. Proceed to download')
					a <- gnomad_fetch_auto(map, type, path, liftover)
				}
			}
			b <- clinvar_fetch(gns, path, liftover)
			progressr::with_progress(c <- lovd3_fetch(gns, path, liftover))
			progressr::with_progress(out <- unique(rbind(a, b, c)))
		}
	)

	return(out)
}

#' Fetch manually downloaded gnomAD variants.
#'
#' Collect known variants for specific genes downloaded manually from gnomAD
#' 
#' @param map Gene data as in map. 
#' @param path Where data were downloaded
#' @param liftover Should variants be lifted over?
#' @param af Maintain allele frequencies?
#' @return vcf data.frame
gnomad_fetch_manual <- function(map, path, liftover, af = FALSE, sig = FALSE){
	if(isTRUE(af) & isTRUE(sig)) stop('Only one of `af`, `sig` should be TRUE')
	pathTmp <- file.path(path, 'gnomad/manual')
	if(!dir.exists(pathTmp)){
		stop('gnomAD files are missing')
		} else {
		message('Fetching gnomad variants')
		gns_tmp <- do.call('rbind', map)$gene_id
		files <- list.files(pathTmp, pattern = paste('gnomAD_v4', gns_tmp, sep = '.*', collapse = '|'), 
			full.names = TRUE)
		if(length(grep(paste(gns_tmp, collapse = '|'), files)) == length(gns_tmp)){
			gnomad_vcf <- do.call('rbind', lapply(as.list(files), function(x){
					tmp <- utils::read.delim(x, sep = ',')
					colnames(tmp)[c(2,3,5,6)] <- c('CHR', 'POS', 'REF', 'ALT')
					tmp$CHR <- paste0('chr', tmp$CHR)
					return(tmp)
			}))
			if(!isTRUE(af) & !isTRUE(sig)){
				gnomad_vcf <- gnomad_vcf[,c('CHR', 'POS', 'REF', 'ALT')]
				} else if (isTRUE(af)){
					gnomad_vcf <- gnomad_vcf[,c('CHR', 'POS', 'Allele.Frequency', 'REF', 'ALT')]
					colnames(gnomad_vcf)[3] <- 'ID'
					gnomad_vcf$ID <- gnomad_vcf$ID*100
					gnomad_vcf <- gnomad_vcf %>% dplyr::mutate(
						ID = dplyr::case_when(
							ID < 0.0001  ~ "<0.0001",
							ID > 0.0001 & ID < 0.001  ~ "<0.001",
							ID > 0.001 & ID < 0.01  ~ "<0.01",
							ID > 0.01 & ID < 0.1  ~ "<0.1",
							.default = as.character(round(ID, 2))
						)  
					)
				} else {
					gnomad_vcf <- gnomad_vcf[,c('CHR', 'POS', 'ClinVar.Clinical.Significance', 'REF', 'ALT')]
					colnames(gnomad_vcf)[3] <- 'ID'
					gnomad_vcf <- na.omit(gnomad_vcf[which(gnomad_vcf$ID %!in% c('', 'not provided')),])
					gnomad_vcf <- gnomad_vcf %>%
						dplyr::mutate(
						  	ID = dplyr::recode(ID, 
								'Conflicting interpretations of pathogenicity' = 'Conflicting',
								'Benign/Likely benign' = 'BLB',
								'Likely benign' = 'LB',
								'Benign' = 'B',
								'Uncertain significance' = 'VUS',
								'Likely pathogenic' = 'LP',
								'Pathogenic/Likely pathogenic' = 'PLP',
								'Pathogenic' = 'P',
								'risk factor' = 'Risk_factor',
								'Conflicting interpretations of pathogenicity; other; risk factor' = 'Conflicting_Risk_factor',
								'drug response' = 'Drug_response',
								'Benign/Likely benign; other' = 'BLB_Other',
								'confers sensitivity' = 'Confers_sensitivity',
								'Benign; confers sensitivity' = 'B_Confers_sensitivity',
								'confers sensitivity; other' = 'Confers_sensitivity_Other'
						  	)
						)
					# To cover other possible entries
					gnomad_vcf$ID <- gsub(' *', '_', gnomad_vcf$ID)
				}
			if(isTRUE(liftover)) gnomad_vcf$POS <- as.data.frame(lift(gnomad_vcf))$start
			return(gnomad_vcf)
		} else {
			stop('gnomAD files are missing')
		}
	}
}

#' Fetch gnomAD variants.
#'
#' Collect known variants for specific genes from gnomAD
#' 
#' @param map Gene data as in map. 
#' @param type Type of gnomAD data to download.
#' @param path Where data were downloaded
#' @param liftover Should variants be lifted over?
#' @param af Maintain allele frequencies?
#' @return vcf data.frame
gnomad_fetch_auto <- function(map, type = c('exomes', 'genomes', 'both'), path, liftover, af = FALSE){
	pathTmp <- file.path(path, 'gnomad')
	suppressWarnings(dir.create(pathTmp, recursive = TRUE))
	type <- match.arg(type)
	message('Fetching gnomad variants')
	gnomad_vcf <- do.call('rbind', lapply(map, function(x){
		switch(type,
			exomes = {
				download_gnomad(paste0(unique(x$seqnames)), 'exomes', pathTmp)
				message('Subsetting gnomAD exomes')
				file.gz <- file.path(pathTmp, paste0('gnomad.exomes.v4.0.sites.', 
					unique(x$seqnames), '.vcf.bgz'))
				file.gz.tbi <- paste(file.gz, ".tbi", sep="")
				gr <- GenomicRanges::GRanges(x$seqnames, IRanges::IRanges(x$start, x$end))
				params <- VariantAnnotation::ScanVcfParam(which = gr)
				vcf <- VariantAnnotation::readVcf(Rsamtools::TabixFile(file.gz), 'hg38', params)
				rg <- SummarizedExperiment::rowRanges(vcf)
				if(!isTRUE(af)) {
					out <- data.frame(
						CHR = as.character(GenomicRanges::seqnames(rg)), 
						POS = as.data.frame(rg@ranges)$start, 
						REF = as.character(S4Vectors::DataFrame(rg)$REF), 
						ALT = as.character(unlist(S4Vectors::DataFrame(rg)$ALT))
					)
					} else {
						inf <- VariantAnnotation::info(vcf)
						out <- data.frame(
							CHR = as.character(GenomicRanges::seqnames(rg)), 
							POS = as.data.frame(rg@ranges)$start, 
							ID = as.character(inf$AF),
							REF = as.character(S4Vectors::DataFrame(rg)$REF), 
							ALT = as.character(unlist(S4Vectors::DataFrame(rg)$ALT))
						)
					}
				if(isTRUE(liftover)) out$POS <- as.data.frame(lift(out))$start
				return(out)
			},
			genomes = {
				download_gnomad(unique(x$seqnames), 'genomes', pathTmp)
				message('Subsetting gnomAD genomes')
				file.gz <- file.path(pathTmp, paste0('gnomad.genomes.v4.0.sites.', 
					unique(x$seqnames), '.vcf.bgz'))
				file.gz.tbi <- paste(file.gz, ".tbi", sep="")
				gr <- GenomicRanges::GRanges(x$seqnames, IRanges::IRanges(x$start, x$end))
				params <- VariantAnnotation::ScanVcfParam(which = gr)
				vcf <- VariantAnnotation::readVcf(Rsamtools::TabixFile(file.gz), 'hg38', params)
				rg <- SummarizedExperiment::rowRanges(vcf)
				if(!isTRUE(af)) {
					out <- data.frame(
						CHR = as.character(GenomicRanges::seqnames(rg)), 
						POS = as.data.frame(rg@ranges)$start, 
						REF = as.character(S4Vectors::DataFrame(rg)$REF), 
						ALT = as.character(unlist(S4Vectors::DataFrame(rg)$ALT))
					)
					} else {
						inf <- VariantAnnotation::info(vcf)
						out <- data.frame(
							CHR = as.character(GenomicRanges::seqnames(rg)), 
							POS = as.data.frame(rg@ranges)$start, 
							ID = as.character(inf$AF),
							REF = as.character(S4Vectors::DataFrame(rg)$REF), 
							ALT = as.character(unlist(S4Vectors::DataFrame(rg)$ALT))
						)
					}
				if(isTRUE(liftover)) out$POS <- as.data.frame(lift(out))$start
				return(out)
			},
			both = {
				download_gnomad(unique(x$seqnames), type, pathTmp)
				out <- list()
				for(k in c('exomes', 'genomes')) {
					message('Subsetting gnomAD ', k)
					file.gz <- file.path(pathTmp, paste0('gnomad.', k, '.v4.0.sites.', 
						unique(x$seqnames), '.vcf.bgz'))
					file.gz.tbi <- paste(file.gz, ".tbi", sep="")
					gr <- GenomicRanges::GRanges(x$seqnames, IRanges::IRanges(x$start, x$end))
					params <- VariantAnnotation::ScanVcfParam(which = gr)
					vcf <- VariantAnnotation::readVcf(Rsamtools::TabixFile(file.gz), 'hg38', params)
					rg <- SummarizedExperiment::rowRanges(vcf)
					if(!isTRUE(af)) {
						out[[k]] <- data.frame(
							CHR = as.character(GenomicRanges::seqnames(rg)), 
							POS = as.data.frame(rg@ranges)$start, 
							REF = as.character(S4Vectors::DataFrame(rg)$REF), 
							ALT = as.character(unlist(S4Vectors::DataFrame(rg)$ALT))
						)
						} else {
							inf <- VariantAnnotation::info(vcf)
							out[[k]] <- data.frame(
								CHR = as.character(GenomicRanges::seqnames(rg)), 
								POS = as.data.frame(rg@ranges)$start, 
								ID = as.character(inf$AF),
								REF = as.character(S4Vectors::DataFrame(rg)$REF), 
								ALT = as.character(unlist(S4Vectors::DataFrame(rg)$ALT))
							)
						}
					if(isTRUE(liftover)) out[[k]]$POS <- as.data.frame(lift(out[[k]]))$start
				}
				return(do.call('rbind', out))
			}

		)
	}))
	rownames(gnomad_vcf) <- NULL
	unlink(c(file.path(pathTmp, 'tmp_gen.vcf'), file.path(pathTmp, 'tmp_ex.vcf')))
	return(gnomad_vcf)
}

#' Fetch ClinVar variants.
#'
#' Collect known variants for specific genes from ClinVar
#' 
#' @param gns Gene names character vector.
#' @param path Where data were downloaded
#' @param liftover Should variants be lifted over?
#' @return vcf data.frame
clinvar_fetch <- function(gns, path, liftover){
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
	# Genes of interest & at least one star
	vcf_sel <- vcf[grep(paste0('^', gns, collapse = '|'), VariantAnnotation::info(vcf)$GENEINFO),]
	rg <- SummarizedExperiment::rowRanges(vcf_sel)
	inf <- VariantAnnotation::info(vcf_sel)
	clinvar_vcf <- data.frame(
		CHR = paste0('chr', as.character(GenomicRanges::seqnames(rg))),
		POS = as.data.frame(rg@ranges)$start, 
		REF = as.character(S4Vectors::DataFrame(rg)$REF), 
		ALT = as.character(unlist(S4Vectors::DataFrame(rg)$ALT))#,
		# GENEINFO = as.character(inf$GENEINFO)
	)
	# system(paste0('bcftools query -f "%CHROM %POS %REF %ALT %GENEINFO\n" ', file.path(pathTmp, 'clinvar.vcf.gz'), 
	# 	' > ', file.path(pathTmp, 'clinvar_tmp.vcf')))
	# clinvar_vcf <- read.table(file.path(pathTmp, 'clinvar_tmp.vcf'))
	# colnames(clinvar_vcf) <- c('CHR', 'POS', 'REF', 'ALT', 'GENEINFO')
	# clinvar_vcf <- clinvar_vcf[grep(paste(gns, collapse = '|'), clinvar_vcf$GENEINFO), 1:4]
	if(isTRUE(liftover)) clinvar_vcf$POS <- as.data.frame(lift(clinvar_vcf))$start
	unlink(file.path(pathTmp, 'clinvar_tmp.vcf'))
	return(clinvar_vcf)
}

#' Fetch lovd3 variants.
#'
#' Collect known variants for specific genes from lovd3
#' 
#' @param gns Gene names character vector.
#' @param path Where data were downloaded
#' @param liftover Should variants be lifted over?
#' @return vcf data.frame
lovd3_fetch <- function(gns, path, liftover){
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
			# message('Now processing ', gene)
			# Genomic data
			start <- grep('\\#\\# Variants_On_Genome', tmp)
			stop <- grep('\\#\\# Variants_On_Transcripts', tmp) -4
			df_genom <- utils::read.delim(i, nrows = stop - (start + 2), skip = start + 3, header = FALSE)
			colnames(df_genom) <- gsub('\\"\\{\\{|\\}\\}\\"', '', strsplit(tmp[start + 3], split = '\t')[[1]])

			# Transcript data
			start <- grep('\\#\\# Variants_On_Transcripts', tmp)
			stop <- grep('\\#\\# Screenings_To_Variants', tmp) -4
			df_trans <- utils::read.delim(i, nrows = stop - (start + 3), skip = start + 4, header = FALSE)
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

		# Remove entries w/o genomic location
		lovd3_vcf <- lovd3_vcf[which(lovd3_vcf$VariantOnGenome.DNA.hg38 %!in% c('', '- ')),]

		# Maintain germline mutations
		lovd3_vcf <- lovd3_vcf[grep('Germline', lovd3_vcf$VariantOnGenome.Genetic_origin),]
		 
		# Remove large duplications and deletions
		lovd3_vcf <- lovd3_vcf[grep('\\)dup|\\)del', lovd3_vcf$VariantOnGenome.DNA.hg38, invert = TRUE),]

		# Remove cases with []
		lovd3_vcf <- lovd3_vcf[grep('\\[',lovd3_vcf$VariantOnGenome.DNA.hg38, invert = TRUE),]

		# Genomic positions
		lovd3_vcf$end <- lovd3_vcf$start <- gsub('g\\.|[a-zA-Z]*>.*$|del.*|dup.*|ins.*|inv.*', '', 
			lovd3_vcf$VariantOnGenome.DNA.hg38)
		lovd3_vcf$start <- suppressWarnings(as.numeric(gsub('_.*$', '', lovd3_vcf$start)))
		lovd3_vcf$end <- suppressWarnings(as.numeric(gsub('^.*_', '', lovd3_vcf$end)))
		lovd3_vcf <- na.omit(lovd3_vcf)

		# Correct some lovd3_vcf typos
		lovd3_vcf$start <- as.numeric(gsub('^0\\.', '', lovd3_vcf$start))
		lovd3_vcf$end <- as.numeric(gsub('^0\\.', '', lovd3_vcf$end))

		# Isolate Ref and Alt bases
		# old2 <- lovd3_vcf
		refAlt <- list()
		message('Processing LOVD3 variants')
		p <- progressr::progressor(steps = nrow(lovd3_vcf))
		for(i in 1:nrow(lovd3_vcf)){
			p(message = sprintf("Adding %g", i))
			# message('Processing ', i, ' of ', nrow(lovd3_vcf))
			tmp <- lovd3_vcf[i,]
			# Correct potential typo
			if(length(grep('^\\.', tmp$VariantOnGenome.DNA.hg38))){
				tmp$VariantOnGenome.DNA.hg38 <- lovd3_vcf[i, 'VariantOnGenome.DNA.hg38'] <- 
					gsub('\\.', 'g\\.', tmp$VariantOnGenome.DNA.hg38)
			}

			if(length(grep('>', tmp$VariantOnGenome.DNA.hg38))){ # Retrieve SNPs

				out <- unlist(strsplit(
					gsub('g\\.[0-9]*', '', tmp$VariantOnGenome.DNA.hg38), split = '>'))
				refAlt[[i]] <- data.frame(REF = out[1], ALT = out[2])

			} else if(length(grep('del$', tmp$VariantOnGenome.DNA.hg38))){ # Retrieve any variant recorded as e.g. `delGGGTGAA`

				# Change coordinates to the previous base
				lovd3_vcf$start[i] <- lovd3_vcf$start[i] - 1
				gseq <- as.character(
					Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, tmp$chromosome, 
						tmp$start - 1, tmp$end)
				)
				refAlt[[i]] <- data.frame(REF = gseq, ALT = unlist(strsplit(gseq, split = ''))[1])

			} else if(length(grep('ins', tmp$VariantOnGenome.DNA.hg38))){ # Gather any del, ins, delins

				# Is the variant recorded as e.g. `delGAinsT`?
				out <- stringr::str_extract(tmp$VariantOnGenome.Published_as, 'del[a-zA-Z]+ins[a-zA-Z]+')

				# If the variant is indeed recorded as e.g. `delGAinsT`
				if(!is.na(out)){
					aa <- unlist(strsplit(out, 'ins'))
					aa[1] <- gsub('del', '', aa[1])
					y <- unlist(strsplit(paste(aa, collapse = ''), ''))
					ifelse(any(y %!in% c('A', 'T', 'C', 'G')),
						out <- data.frame(REF = NA, ALT = NA),
						out <- data.frame(REF = aa[1], ALT = aa[2])
					)
				}

				# Is the variant recorded as e.g. `insT`?
				if(length(out) == 1 && is.na(out)){
					# Change coordinates to the previous base
					lovd3_vcf$start[i] <- lovd3_vcf$start[i] - 1

					# Isolate insertion
					out <- stringr::str_extract(tmp$VariantOnGenome.Published_as, 'ins[a-zA-Z]+')
					out <- gsub('ins', '', out)

					# Previous base
					ref <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, tmp$chromosome, 
						tmp$start - 1, tmp$start - 1)

					y <- unlist(strsplit(out, ''))
					ifelse(any(y %!in% c('A', 'T', 'C', 'G')),
						out <- data.frame(REF = NA, ALT = NA),
						out <- data.frame(REF = ref, ALT = paste0(ref, out))
					)

				}

				# Is the variant recorded as e.g. `delinsTGC` in `VariantOnGenome.DNA.hg38`?
				if(length(out) == 1 && is.na(out)){
					out <- stringr::str_extract(tmp$VariantOnGenome.DNA.hg38, 'delins[a-zA-Z]+')
					out <- gsub('delins', '', out)
					y <- unlist(strsplit(out, ''))

					# Previous base
					ref <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, tmp$chromosome, 
						tmp$start, tmp$end)

					ifelse(any(y %!in% c('A', 'T', 'C', 'G')),
						out <- data.frame(REF = NA, ALT = NA),
						out <- data.frame(REF = ref, ALT = out)
					)
				}

				# Is the variant recorded as e.g. `insTGC` in `VariantOnGenome.DNA.hg38`?
				if(length(out) == 1 && is.na(out)){

					# Change coordinates to the previous base
					lovd3_vcf$start[i] <- lovd3_vcf$start[i] - 1

					# Isolate insertion
					out <- stringr::str_extract(tmp$VariantOnGenome.DNA.hg38, 'ins[a-zA-Z]+')
					out <- gsub('ins', '', out)

					# Previous base
					ref <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, tmp$chromosome, 
						tmp$start - 1, tmp$start - 1)

					y <- unlist(strsplit(out, ''))
					ifelse(any(y %!in% c('A', 'T', 'C', 'G')),
						out <- data.frame(REF = NA, ALT = NA),
						out <- data.frame(REF = ref, ALT = paste0(ref, out))
					)
				}

				# If all the above fail
				if(length(out) == 1 && is.na(out)){
					out <- data.frame(REF = NA, ALT = NA)
				}

				refAlt[[i]] <- out

			} else if(length(grep('dup', tmp$VariantOnGenome.DNA.hg38))){ # Gather any dups

				out <- stringr::str_extract(tmp$VariantOnGenome.Published_as, 'dup[a-zA-Z]+')

				# If the variant is indeed recorded as e.g. `478_480dupGGA`
				if(!is.na(out)){
					aa <- gsub('dup', '', out)
					y <- unlist(strsplit(paste(aa, collapse = ''), ''))
					ifelse(any(y %!in% c('A', 'T', 'C', 'G')),
						out <- data.frame(REF = NA, ALT = NA),
						out <- data.frame(REF = aa, ALT = paste(rep(aa, 2), collapse = ''))
					)
				}

				# Is the variant recorded as e.g. `1280_1281insA`?
				if(length(out) == 1 && is.na(out)){
					
					# Change coordinates to the previous base
					lovd3_vcf$start[i] <- lovd3_vcf$start[i] - 1

					# Isolate insertion
					out <- stringr::str_extract(tmp$VariantOnGenome.Published_as, 'ins[a-zA-Z]+')
					out <- gsub('ins', '', out)

					# Previous base
					ref <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, tmp$chromosome, 
						tmp$start - 1, tmp$start - 1)

					y <- unlist(strsplit(out, ''))
					ifelse(any(y %!in% c('A', 'T', 'C', 'G')),
						out <- data.frame(REF = NA, ALT = NA),
						out <- data.frame(REF = ref, ALT = paste0(ref, out))
					)

				}

				# If all the above fail
				if(length(out) == 1 && is.na(out)){
					out <- data.frame(REF = NA, ALT = NA)
				}

				refAlt[[i]] <- out

			} else { # A handfull of exceptions
				refAlt[[i]] <- data.frame(REF = NA, ALT = NA)
			}
		}
		lovd3_vcf <- cbind(lovd3_vcf, do.call('rbind', refAlt))

		# Keep important columns
		lovd3_vcf <- unique(lovd3_vcf[,c('chromosome', 'start', 'REF', 'ALT')])
		rownames(lovd3_vcf) <- NULL
		lovd3_vcf <- na.omit(lovd3_vcf)
		lovd3_vcf <- lovd3_vcf[grep('Y|W|K|D', lovd3_vcf$ALT, invert = TRUE),]
		lovd3_vcf <- lovd3_vcf[grep('Y|W|K|D', lovd3_vcf$ALT, invert = TRUE),]
		colnames(lovd3_vcf)[1:2] <- c('CHR', 'POS')

		} else {
			lovd3_vcf <- data.frame(CHR = character(), POS = character(),
				REF = character(), ALT = character())
		}
	if(isTRUE(liftover)) lovd3_vcf$POS <- as.data.frame(lift(lovd3_vcf))$start
	unlink(file.path(pathTmp, '*'))
	return(lovd3_vcf)
}
