#' Create annotation files for IonReporter.
#'
#' A wrapper function for downloading and organizing annotation files.
#'
#' @param gns Gene names character vector.
#' @param annotators Annotators to create files for.
#' @param panelName Path and prefix for output file.  
#' @param path Where database files will be downloaded.  
#' @param type Type of gnomAD data to download. Applicable only for `intervar` and `all` annotators
#' @return IonReporter annotation files.
#' @export
annotate <- function(gns, annotators = c('revel', 'alphamissense', 'intervar', 'all'),
	panelName = './results/panel', path = './dbs', type = c('exomes', 'genomes', 'both')){

	annotators <- match.arg(annotators)
	suppressWarnings(dir.create('./results'))
	switch(annotators,
		revel = {
			vcf_body <- constructREVEL(gns, path)
			vcf_header <- c(
				"##fileformat=VCFv4.2",
				"##HITLEVEL=allele",
				"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample"
			)
			write.table(vcf_header, file = paste0(panelName, '_revel.vcf'), 
				sep = '\t', quote = FALSE, col.names = FALSE, 
				row.names = FALSE)
			write.table(vcf_body, file =  paste0(panelName, '_revel.vcf'), 
				append = TRUE, sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE, 
				na = '')
		},
		alphamissense = {
			vcf_body <- constructAM(gns, path)
			vcf_header <- c(
				"##fileformat=VCFv4.2",
				"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
				"##HITLEVEL=genotype",
				"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample"
			)
			write.table(vcf_header, file = paste0(panelName, '_alphaMissense.vcf'), 
				sep = '\t', quote = FALSE, col.names = FALSE, 
				row.names = FALSE)
			write.table(vcf_body, file = paste0(panelName, '_alphaMissense.vcf'), 
				append = TRUE, sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE, 
				na = '')
		},
		intervar = {
			message('Collecting variants')
			vars <- collectVars(gns, databases = 'all', path = path, type = type)
			message('Pathogenicity prediction')
			progressr::with_progress(vars <- annotateInterVar(vars))
			vars[is.na(vars)] <- NULL
			message(length(vars), ' variants were annotated.')
			vcf_body <- do.call('rbind', lapply(vars, function(x){
				tmp <- as.data.frame(x)
				crit <- tmp[,8:ncol(tmp)]
				data.frame(CHR = paste0('chr', tmp$Chromosome),
					POS = tmp$Position, 
					INFO = paste0(tmp$Intervar, ' (', 
						paste(colnames(crit)[crit == 1], collapse = '; '),
						')'),
					REF = tmp$Ref_allele, ALT = tmp$Alt_allele, QUAL = '.',
					FILTER = 'PASS', INFO = '.', FORMAT = '.', Sample = '.')
				}))
			vcf_header <- c(
				"##fileformat=VCFv4.2",
				"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
				"##HITLEVEL=genotype",
				"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample"
			)
			write.table(vcf_header, file = paste0(panelName, '_intervar.vcf'), sep = '\t', quote = FALSE, col.names = FALSE, 
				row.names = FALSE)
			write.table(vcf_body, file = paste0(panelName, '_intervar.vcf'), append = TRUE, sep = '\t', quote = FALSE, 
				col.names = FALSE, row.names = FALSE, na = '')
		},
		all = {
			# ===================================
			# REVEL
			# ===================================
			vcf_body <- constructREVEL(gns, path)
			vcf_header <- c(
				"##fileformat=VCFv4.2",
				"##HITLEVEL=allele",
				"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample"
			)
			write.table(vcf_header, file = paste0(panelName, '_revel.vcf'), 
				sep = '\t', quote = FALSE, col.names = FALSE, 
				row.names = FALSE)
			write.table(vcf_body, file =  paste0(panelName, '_revel.vcf'), 
				append = TRUE, sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE, 
				na = '')

			# ===================================
			# AlphaMissense
			# ===================================
			vcf_body <- constructAM(gns, path)
			vcf_header <- c(
				"##fileformat=VCFv4.2",
				"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
				"##HITLEVEL=genotype",
				"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample"
			)
			write.table(vcf_header, file = paste0(panelName, '_alphaMissense.vcf'), 
				sep = '\t', quote = FALSE, col.names = FALSE, 
				row.names = FALSE)
			write.table(vcf_body, file = paste0(panelName, '_alphaMissense.vcf'), 
				append = TRUE, sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE, 
				na = '')

			# ===================================
			# InterVar
			# ===================================
			vars <- collectVars(gns, databases = 'all', path)
			progressr::with_progress(vars <- annotateInterVar(vars))
			vars[is.na(vars)] <- NULL
			message(length(vars), ' variants were annotated.')
			vcf_body <- do.call('rbind', lapply(vars, function(x){
				tmp <- as.data.frame(x)
				crit <- tmp[,8:ncol(tmp)]
				data.frame(CHR = paste0('chr', tmp$Chromosome),
					POS = tmp$Position, 
					INFO = paste0(tmp$Intervar, ' (', 
						paste(colnames(crit)[crit == 1], collapse = '; '),
						')'),
					REF = tmp$Ref_allele, ALT = tmp$Alt_allele, QUAL = '.',
					FILTER = 'PASS', INFO = '.', FORMAT = '.', Sample = '.')
				}))
			vcf_header <- c(
				"##fileformat=VCFv4.2",
				"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
				"##HITLEVEL=genotype",
				"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample"
			)
			write.table(vcf_header, file = paste0(panelName, '_intervar.vcf'), sep = '\t', quote = FALSE, col.names = FALSE, 
				row.names = FALSE)
			write.table(vcf_body, file = paste0(panelName, '_intervar.vcf'), append = TRUE, sep = '\t', quote = FALSE, 
				col.names = FALSE, row.names = FALSE, na = '')
		}
	)
}