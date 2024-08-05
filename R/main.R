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
#' @param liftover Should variants be lifted over?
#' @param saveRaw Save all variants without any annotation. Applicable only for `intervar` and `all` annotators.
#' @return IonReporter annotation files.
#' @export
annotate <- function(gns, annotators = c('revel', 'alphamissense', 'clinvar_sig', 'gnomad_af', 'intervar', 'all'),
	databases = c('gnomad_man', 'gnomad_auto', 'clinvar', 'lovd3', 'all'), panelName = './results/panel',
	path = './dbs', type = c('exomes', 'genomes', 'both'), liftover = FALSE, saveRaw = FALSE){

	suppressWarnings(dir.create('./results'))
	annotators <- match.arg(annotators)
	databases <- match.arg(databases)
	type <- match.arg(type)
	switch(annotators,
		revel = {
			message('\n=========================\nREVEL annotation\n=========================')
			vcf_body <- constructREVEL(gns, path, liftover)
			utils::write.table(vcf_header_allele, file = paste0(panelName, '_revel.vcf'), 
				sep = '\t', quote = FALSE, col.names = FALSE, 
				row.names = FALSE)
			utils::write.table(vcf_body, file =  paste0(panelName, '_revel.vcf'), 
				append = TRUE, sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE, 
				na = '')
		},
		alphamissense = {
			message('\n=========================\nAlphaMissense annotation\n=========================')
			vcf_body <- constructAM(gns, path, liftover)
			utils::write.table(vcf_header_allele, file = paste0(panelName, '_alphaMissense.vcf'), 
				sep = '\t', quote = FALSE, col.names = FALSE, 
				row.names = FALSE)
			utils::write.table(vcf_body, file = paste0(panelName, '_alphaMissense.vcf'), 
				append = TRUE, sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE, 
				na = '')
		},
		clinvar_sig = {
			message('\n========================\nClinVar annotation\n========================')
			message('ClinVar annotation is currently retrieved from gnomAD')
			vcf_body <- constructClinVar(gns, path, liftover)
			utils::write.table(vcf_header_allele, file = paste0(panelName, '_clinvar_sig.vcf'), 
				sep = '\t', quote = FALSE, col.names = FALSE, 
				row.names = FALSE)
			utils::write.table(vcf_body, file =  paste0(panelName, '_clinvar_sig.vcf'), 
				append = TRUE, sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE, 
				na = '')
		},
		gnomad_af = {
			message('\n========================\ngnomAD frequencies\n========================')
			vcf_body <- constructAF(gns, path, databases, type, liftover)
			utils::write.table(vcf_header_allele, file = paste0(panelName, '_gnomad_af.vcf'), 
				sep = '\t', quote = FALSE, col.names = FALSE, 
				row.names = FALSE)
			utils::write.table(vcf_body, file =  paste0(panelName, '_gnomad_af.vcf'), 
				append = TRUE, sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE, 
				na = '')
		},
		intervar = {
			message('\n========================\nInterVar annotation\n========================')
			vars <- collectVars(gns, databases, path, type, liftover)
			if(isTRUE(saveRaw)){
				varsOut <- data.frame(vars[,1:2], ID = '.', vars[,3:4], QUAL = '.', FILTER = 'PASS',
					INFO = '.', FORMAT = '.', Sample = '.')
				utils::write.table(varsOut, file = paste0(panelName, '_raw.txt'), sep = '\t', 
					quote = FALSE, row.names = FALSE)
			}
			if(nrow(vars) != 0){
				message('\n========================\nPathogenicity prediction\n========================')
				progressr::with_progress(vars <- annotateInterVar(vars, liftover))
				vars[is.na(vars)] <- NULL
				message(length(vars), ' variants were successfully annotated.')
				vcf_body_hom <- vcf_body_het <- do.call('rbind', lapply(vars, function(x){
					tmp <- as.data.frame(x)
					crit <- tmp[,8:ncol(tmp)]
					data.frame(CHR = paste0('chr', tmp$Chromosome),
						POS = tmp$Position, 
						ID = paste0(gsub(' ', '_', tmp$Intervar), '(', 
							paste(colnames(crit)[crit == 1], collapse = ';'),
							')'),
						REF = tmp$Ref_allele, ALT = tmp$Alt_allele, QUAL = '.',
						FILTER = 'PASS', INFO = '.', FORMAT = 'GT:GQ:DP:AD', Sample = '.')
					
				}))
				gq_dp <- sample(80:300, 3)
				vcf_body_het$Sample <- paste('0/1', gq_dp[1], gq_dp[2]+gq_dp[3], paste(gq_dp[2], gq_dp[3], sep = ','), sep = ':')
				vcf_body_hom$Sample <- paste('1/1', gq_dp[1], gq_dp[2]+gq_dp[3], paste(gq_dp[2], gq_dp[3], sep = ','), sep = ':')
				} else {
					message('Pathogenicity prediction aborted.\nNo variants found.')
					vcf_body_het <- vcf_body_hom <- data.frame(CHR = character(), POS = character(),
						REF = character(), ALT = character())
				}
			for(i in c('het', 'hom')){
				utils::write.table(vcf_header_genotype, file = paste(panelName, databases, i, 'intervar.vcf', sep = '_'), 
					sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)
				utils::write.table(get(paste0('vcf_body_', i)), file = paste(panelName, databases, i, 'intervar.vcf', sep = '_'), append = TRUE, sep = '\t', quote = FALSE, 
					col.names = FALSE, row.names = FALSE, na = '')
			}
		},
		all = {
			message('\n=========================\nREVEL annotation\n=========================')
			vcf_body <- constructREVEL(gns, path, liftover)
			utils::write.table(vcf_header_allele, file = paste0(panelName, '_revel.vcf'), 
				sep = '\t', quote = FALSE, col.names = FALSE, 
				row.names = FALSE)
			utils::write.table(vcf_body, file =  paste0(panelName, '_revel.vcf'), 
				append = TRUE, sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE, 
				na = '')

			message('\n=========================\nAlphaMissense annotation\n=========================')
			vcf_body <- constructAM(gns, path, liftover)
			utils::write.table(vcf_header_allele, file = paste0(panelName, '_alphaMissense.vcf'), 
				sep = '\t', quote = FALSE, col.names = FALSE, 
				row.names = FALSE)
			utils::write.table(vcf_body, file = paste0(panelName, '_alphaMissense.vcf'), 
				append = TRUE, sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE, 
				na = '')

			message('\n========================\nClinVar annotation\n========================')
			vcf_body <- constructClinVar(gns, path, liftover)
			utils::write.table(vcf_header_allele, file = paste0(panelName, '_clinvar_sig.vcf'), 
				sep = '\t', quote = FALSE, col.names = FALSE, 
				row.names = FALSE)
			utils::write.table(vcf_body, file =  paste0(panelName, '_clinvar_sig.vcf'), 
				append = TRUE, sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE, 
				na = '')

			message('\n========================\ngnomAD frequencies\n========================')
			vcf_body <- constructAF(gns, path, databases, type, liftover)
			utils::write.table(vcf_header_allele, file = paste0(panelName, '_gnomad_af.vcf'), 
				sep = '\t', quote = FALSE, col.names = FALSE, 
				row.names = FALSE)
			utils::write.table(vcf_body, file =  paste0(panelName, '_gnomad_af.vcf'), 
				append = TRUE, sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE, 
				na = '')

			message('\n========================\nInterVar annotation\n========================')
			vars <- collectVars(gns, databases, path, type, liftover)
			if(isTRUE(saveRaw)){
				varsOut <- data.frame(vars[,1:2], ID = '.', vars[,3:4], QUAL = '.', FILTER = 'PASS',
					INFO = '.', FORMAT = '.', Sample = '.')
				utils::write.table(varsOut, file = paste0(panelName, '_raw.txt'), sep = '\t', 
					quote = FALSE, row.names = FALSE)
			}
			if(nrow(vars) != 0){
				message('\n========================\nPathogenicity prediction\n========================')
				progressr::with_progress(vars <- annotateInterVar(vars, liftover))
				vars[is.na(vars)] <- NULL
				message(length(vars), ' variants were successfully annotated.')
				vcf_body_hom <- vcf_body_het <- do.call('rbind', lapply(vars, function(x){
					tmp <- as.data.frame(x)
					crit <- tmp[,8:ncol(tmp)]
					data.frame(CHR = paste0('chr', tmp$Chromosome),
						POS = tmp$Position, 
						ID = paste0(gsub(' ', '_', tmp$Intervar), '(', 
							paste(colnames(crit)[crit == 1], collapse = ';'),
							')'),
						REF = tmp$Ref_allele, ALT = tmp$Alt_allele, QUAL = '.',
						FILTER = 'PASS', INFO = '.', FORMAT = 'GT:GQ:DP:AD', Sample = '.')
				}))
				gq_dp <- sample(80:300, 3)
				vcf_body_het$Sample <- paste('0/1', gq_dp[1], gq_dp[2]+gq_dp[3], paste(gq_dp[2], gq_dp[3], sep = ','), sep = ':')
				vcf_body_hom$Sample <- paste('1/1', gq_dp[1], gq_dp[2]+gq_dp[3], paste(gq_dp[2], gq_dp[3], sep = ','), sep = ':')
				} else {
					message('Pathogenicity prediction aborted.\nNo variants found.')
					vcf_body_het <- vcf_body_hom <- data.frame(CHR = character(), POS = character(),
						REF = character(), ALT = character())
				}
			for(i in c('het', 'hom')){
				utils::write.table(vcf_header_genotype, file = paste(panelName, databases, i, 'intervar.vcf', sep = '_'), 
					sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)
				utils::write.table(get(paste0('vcf_body_', i)), file = paste(panelName, databases, i, 'intervar.vcf', sep = '_'), append = TRUE, sep = '\t', quote = FALSE, 
					col.names = FALSE, row.names = FALSE, na = '')
			}
		}
	)
}