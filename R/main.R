#' Create annotation files for IonReporter.
#'
#' A wrapper function for downloading and organizing annotation files.
#'
#' @param gns Gene names character vector.
#' @param annotators Annotators to create files for.
#' @param databases Databases to fetch annotation from. 
#' @param panelName Path and prefix for output file.  
#' @param path Where database files will be downloaded.  
#' @param type Type of gnomAD data to download. Applicable only for `intervar` and `all` annotators.
#' @param saveRaw Save all variants without any annotation. Applicable only for `intervar` and `all` annotators.
#' @return IonReporter annotation files.
#' @export
annotate <- function(gns, annotators = c('revel', 'alphamissense', 'clinvar_sig', 'intervar', 'all'),
	databases = c('gnomad_man', 'gnomad_auto', 'clinvar', 'lovd3', 'all'), panelName = './results/panel',
	path = './dbs', type = c('exomes', 'genomes', 'both'), saveRaw = FALSE){

	suppressWarnings(dir.create('./results'))
	annotators <- match.arg(annotators)
	databases <- match.arg(databases)
	type <- match.arg(type)
	switch(annotators,
		revel = {
			vcf_body <- constructREVEL(gns, path)
			utils::write.table(vcf_header_allele, file = paste0(panelName, '_revel.vcf'), 
				sep = '\t', quote = FALSE, col.names = FALSE, 
				row.names = FALSE)
			utils::write.table(vcf_body, file =  paste0(panelName, '_revel.vcf'), 
				append = TRUE, sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE, 
				na = '')
		},
		alphamissense = {
			vcf_body <- constructAM(gns, path)
			utils::write.table(vcf_header_genotype, file = paste0(panelName, '_alphaMissense.vcf'), 
				sep = '\t', quote = FALSE, col.names = FALSE, 
				row.names = FALSE)
			utils::write.table(vcf_body, file = paste0(panelName, '_alphaMissense.vcf'), 
				append = TRUE, sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE, 
				na = '')
		},
		clinvar_sig = {
			vcf_body <- constructClinVar(gns, path)
			utils::write.table(vcf_header_allele, file = paste0(panelName, '_clinvar_sig.vcf'), 
				sep = '\t', quote = FALSE, col.names = FALSE, 
				row.names = FALSE)
			utils::write.table(vcf_body, file =  paste0(panelName, '_clinvar_sig.vcf'), 
				append = TRUE, sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE, 
				na = '')
		},
		intervar = {
			message('Collecting variants')
			vars <- collectVars(gns, databases, path, type)
			if(isTRUE(saveRaw)){
				varsOut <- data.frame(vars[,1:2], ID = '.', vars[,3:4], QUAL = '.', FILTER = 'PASS',
					INFO = '.', FORMAT = '.', Sample = '.')
				utils::write.table(varsOut, file = paste0(panelName, '_raw.txt'), sep = '\t', 
					quote = FALSE, row.names = FALSE)
			}
			if(nrow(vars) != 0){
				message('Pathogenicity prediction')
				progressr::with_progress(vars <- annotateInterVar(vars))
				vars[is.na(vars)] <- NULL
				message(length(vars), ' variants were annotated.')
				vcf_body <- do.call('rbind', lapply(vars, function(x){
					tmp <- as.data.frame(x)
					crit <- tmp[,8:ncol(tmp)]
					data.frame(CHR = paste0('chr', tmp$Chromosome),
						POS = tmp$Position, 
						ID = paste0(tmp$Intervar, ' (', 
							paste(colnames(crit)[crit == 1], collapse = '; '),
							')'),
						REF = tmp$Ref_allele, ALT = tmp$Alt_allele, QUAL = '.',
						FILTER = 'PASS', INFO = '.', FORMAT = '.', Sample = '.')
					}))
				} else {
					message('Pathogenicity prediction aborted.\nNo variants found.')
					vcf_body <- data.frame(CHR = character(), POS = character(),
						REF = character(), ALT = character())
				}
			utils::write.table(vcf_header_genotype, file = paste(panelName, databases, 'intervar.vcf', sep = '_'), 
				sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)
			utils::write.table(vcf_body, file = paste0(panelName, '_intervar.vcf'), append = TRUE, sep = '\t', quote = FALSE, 
				col.names = FALSE, row.names = FALSE, na = '')
		},
		all = {
			# ===================================
			# REVEL
			# ===================================
			vcf_body <- constructREVEL(gns, path)
			utils::write.table(vcf_header_allele, file = paste0(panelName, '_revel.vcf'), 
				sep = '\t', quote = FALSE, col.names = FALSE, 
				row.names = FALSE)
			utils::write.table(vcf_body, file =  paste0(panelName, '_revel.vcf'), 
				append = TRUE, sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE, 
				na = '')

			# ===================================
			# AlphaMissense
			# ===================================
			vcf_body <- constructAM(gns, path)
			utils::write.table(vcf_header_genotype, file = paste0(panelName, '_alphaMissense.vcf'), 
				sep = '\t', quote = FALSE, col.names = FALSE, 
				row.names = FALSE)
			utils::write.table(vcf_body, file = paste0(panelName, '_alphaMissense.vcf'), 
				append = TRUE, sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE, 
				na = '')

			# ===================================
			# InterVar
			# ===================================
			message('Collecting variants')
			vars <- collectVars(gns, databases, path, type)
			if(isTRUE(saveRaw)){
				varsOut <- data.frame(vars[,1:2], ID = '.', vars[,3:4], QUAL = '.', FILTER = 'PASS',
					INFO = '.', FORMAT = '.', Sample = '.')
				utils::write.table(varsOut, file = paste0(panelName, '_raw.txt'), sep = '\t', 
					quote = FALSE, row.names = FALSE)
			}
			if(nrow(vars) != 0){
				message('Pathogenicity prediction')
				progressr::with_progress(vars <- annotateInterVar(vars))
				vars[is.na(vars)] <- NULL
				message(length(vars), ' variants were annotated.')
				vcf_body <- do.call('rbind', lapply(vars, function(x){
					tmp <- as.data.frame(x)
					crit <- tmp[,8:ncol(tmp)]
					data.frame(CHR = paste0('chr', tmp$Chromosome),
						POS = tmp$Position, 
						ID = paste0(tmp$Intervar, ' (', 
							paste(colnames(crit)[crit == 1], collapse = '; '),
							')'),
						REF = tmp$Ref_allele, ALT = tmp$Alt_allele, QUAL = '.',
						FILTER = 'PASS', INFO = '.', FORMAT = '.', Sample = '.')
					}))
				} else {
					message('Pathogenicity prediction aborted.\nNo variants found.')
					vcf_body <- data.frame(CHR = character(), POS = character(),
						REF = character(), ALT = character())
				}
			utils::write.table(vcf_header_genotype, file = paste(panelName, databases, 'intervar.vcf', sep = '_'), 
				sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)
			utils::write.table(vcf_body, file = paste0(panelName, '_intervar.vcf'), append = TRUE, sep = '\t', quote = FALSE, 
				col.names = FALSE, row.names = FALSE, na = '')
		}
	)
}