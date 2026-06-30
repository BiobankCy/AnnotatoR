#' Handling errors
#'
#' @param j json file
#' @return LOVD3 gene files.
process_json <- function(j) {tryCatch({
	rjson::fromJSON(j)
}, error = function(e) {
	cat("ERROR :", conditionMessage(e), "\n")
	NA
})}

#' Annotate via InterVar's API.
#'
#' Predict variants' pathogenicity according to
#' ACMG criteria as implemented by InterVar.
#'
#' @param vcf A variant table as the one created by `collectVars`
#' @param liftover Were variants lifted over?
#' @return A nested list with InterVar's ACMG and other annotation pieces per variant. 
#' @export
annotateInterVar <- function(vcf, liftover){
	p <- progressr::progressor(steps = nrow(vcf))
	build <- ifelse(liftover, 'hg19', 'hg38')
	out <- list()
	for(i in 1:nrow(vcf)){
		p(message = sprintf("Adding %g", i))
		tmp <- suppressMessages(httr::content(httr::GET(paste0(
			'http://wintervar.wglab.org/api_new.php?queryType=position&chr=', gsub('chr', '', vcf[i, 'CHR']), 
			'&pos=', vcf[i, 'POS'],
			'&ref=', vcf[i, 'REF'],
			'&alt=', vcf[i, 'ALT'],
			paste0('&build=', build)
		)), 'text'))

		if(tmp == ''){
			out[[i]] <- NA	
		} else {
		 	out[[i]] <- process_json(tmp) 
		}
	}
	return(out)
}

#' Annotate genes using Ensembl data
#'
#' @param ensdb An EnsDb package.
#' @param gns Gene names (HGNC) character vector.
#' @param trans Should transcript level data be returned?
#' @return A data frame with genes (default) or ensembl transcripts as rows and Ensembl annotation as columns.
#' @export
map_fetch <- function(ensdb, gns, trans = FALSE){
	if(isTRUE(trans)){
		edb <- getExportedValue(ensdb, ensdb)
		ensembldb::seqlevelsStyle(edb) <- "UCSC"
		suppressWarnings(granges <- ensembldb::genes(edb))
		suppressWarnings(tranges <- ensembldb::transcripts(edb))
		gns_id <- unique(as.data.frame(granges[GenomicRanges::elementMetadata(granges)$symbol %in% gns])[,'gene_id'])
		gns_id <- gns_id[grep('ENSG', gns_id)]
		map <- as.data.frame(tranges[GenomicRanges::elementMetadata(tranges)$gene_id %in% gns_id])
		map$seqnames <- as.character(map$seqnames)
		map <- map[which(map$seqnames %in% paste0('chr', c(as.character(1:22), 'M', 'X', 'Y'))), ]
		} else {
			edb <- getExportedValue(ensdb, ensdb)
			ensembldb::seqlevelsStyle(edb) <- "UCSC"
			suppressWarnings(granges <- ensembldb::genes(edb))
			map <- as.data.frame(granges[GenomicRanges::elementMetadata(granges)$symbol %in% gns])
			map <- map[grep('ENSG', map$gene_id), ]
			map$seqnames <- as.character(map$seqnames)
			map <- map[which(map$seqnames %in% paste0('chr', c(as.character(1:22), 'M', 'X', 'Y'))), ]
		}
	return(map)
}

#' LiftOver hg38 to hg19 genomic coordinates using rtracklayer functionalities
#'
#' @param vcf A vcf dataframe with CHR and POS columns
#' @param ... Arguments for methods
#' @return A GRanges object of lifted over genomic coordinates
lift <- function(vcf, ...){
	if('index' %in% names(vcf)){ # maybe redundant
		gr38 <- GenomicRanges::GRanges(paste(vcf$CHR, vcf$POS, sep = ':'), index = vcf$index)
	} else {
		gr38 <- GenomicRanges::GRanges(paste(vcf$CHR, vcf$POS, sep = ':'))
	}
	ch <- rtracklayer::import.chain(system.file(package = 'liftOver', 'extdata', 'hg38ToHg19.over.chain'))
	gr19 <- unlist(rtracklayer::liftOver(gr38, ch))
}

#' Annotate via MutationTaster's API.
#'
#' Retrieve variants' pathogenicity predictions from MutationTaster.
#'
#' @param vcf A variant table as the one created by `collectVars`.
#' @param liftover Whether variants have been lifted over (boolean).
#' @return Annotated variants with MutationTaster2021 and MutationTaster2025 predictions for GRCh37 (liftover = TRUE) and GRCh38 (liftover = FALSE) variants, respectively. 
#' @export
taste <- function(vcf, liftover){
	p <- progressr::progressor(steps = nrow(vcf))
	out <- list()
	mainAddress <- ifelse(isTRUE(liftover), 'https://www.genecascade.org/MT2021/MT_API102.cgi?variants=',
		'https://www.genecascade.org/MutationTaster2025/modperl/API.cgi?variants=')
	for(i in 1:nrow(vcf)){
		p(message = sprintf("Adding %g", i))
		mt <- httr::GET(paste0(mainAddress, vcf[i, 'CHR'], ':', vcf[i, 'POS'], vcf[i, 'REF'], '%3E', vcf[i, 'ALT']))

		# Check status code & ERROR response
		if(mt$status_code == 414 | grepl('^ERROR', rawToChar(mt$content))){
			out[[i]] <- as.data.frame(matrix(NA, nrow = 0, ncol = 15))
			} else {
			cont <- lapply(strsplit(rawToChar(mt$content), split = '\\n'), strsplit, split = '\\t')[[1]]
			hd <- cont[[1]]
			hd[2:5] <- toupper(hd[2:5])
			if(length(cont) > 1){
				out[[i]] <- do.call('rbind', lapply(cont[2:length(cont)], function(x){
					if(length(x) == 14) x <- c(x, '')
					return(x)
					}))
				} else {
					out[[i]] <- as.data.frame(matrix(NA, nrow = 0, ncol = 15))
				}

			}
			colnames(out[[i]]) <- hd
	}
	out <- do.call('rbind', out)
	# Correct chromosome notation
	out$CHR <- paste0('chr', out$CHR)
	out <- out %>%
		dplyr::mutate(
		  	CHR = dplyr::recode(CHR, 
				'chr23' = 'chrX',
				'chr24' = 'chrY'
		  	)
		)
	out <- merge(vcf, out, by = c('CHR', 'POS', 'REF', 'ALT'))
	return(out)
}
