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
#' @return Annotated variants
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
#' @param gns Gene names character vector.
#' @param trans Should return transcript data?
#' @return Annotated genes
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

#' LiftOver hg38 to hg19
#'
#' @param vcf A vcf dataframe with CHR and POS columns
#'
#' @return Lifted over positions
lift <- function(vcf, ...){
	gr38 <- GenomicRanges::GRanges(paste(vcf$CHR, vcf$POS, sep = ':'))
	ch <- rtracklayer::import.chain(system.file(package = 'liftOver', 'extdata', 'hg38ToHg19.over.chain'))
	gr19 <- unlist(rtracklayer::liftOver(gr38, ch))
}
