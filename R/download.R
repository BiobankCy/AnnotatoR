#' Download gnomad chromosome files
#'
#' Downloads gnomAD vcf files according to specified chromosome.
#'
#' @param chr Chromosome name
#' @param type Data types to download
#' @param dir Where to download data
#' @return gnomAD vcf files
#' @export
download_gnomad <- function(chr, type = c('exomes', 'genomes', 'both'), 
	path){
	suppressWarnings(dir.create(path, recursive = TRUE))
	type <- match.arg(type)
	switch(type,
		exomes = {
			f <- paste0('gnomad.exomes.v4.0.sites.', chr, '.vcf.bgz')
			if(file.exists(f)) {
				return
				} else  {
					lnk <- paste0(' https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/exomes/gnomad.exomes.v4.0.sites.', 
						chr, '.vcf.bgz')
					system(paste0("curl -sI", lnk, " | grep -i Content-Length | cut -d ' ' -f 2 | uniq > ", path, "/size.txt"))
					message('You are about to download ~', round(as.integer(readLines(file.path(path, 'size.txt')))*10^-9, 2), 'GiB of data.') 
					if(isTRUE(askYesNo('Do you want to proceed?'))){
						system(paste0('wget -q --show-progress -P ', path, lnk))
						lnkIndex <- paste0(lnk, '.tbi')
						system(paste0('wget -q --show-progress -P ', path, lnkIndex))
						unlink(c(file.path(path, 'size.txt'), 'NUL'))
					} else {
						unlink(c(file.path(path, 'size.txt'), 'NUL'))
						stop('Select another annotator set.')
					}
				}
		},
		genomes = {
			f <- paste0('gnomad.exomes.v4.0.sites.', chr, '.vcf.bgz')
			if(file.exists(f)) {
				return
				} else  {
					lnk <- paste0(' https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/genomes/gnomad.genomes.v4.0.sites.', 
						chr, '.vcf.bgz')	
					system(paste0("curl -sI", lnk, " | grep -i Content-Length | cut -d ' ' -f 2 | uniq > ", path, "/size.txt"))
					message('You are about to download ~', round(as.integer(readLines(file.path(path, 'size.txt')))*10^-9, 2), 'GiB of data.') 
					if(isTRUE(askYesNo('Do you want to proceed?'))){
						system(paste0('wget -q --show-progress -P ', path, lnk))
						lnkIndex <- paste0(lnk, '.tbi')
						system(paste0('wget -q --show-progress -P ', path, lnkIndex))
						unlink(c(file.path(path, 'size.txt'), 'NUL'))
					} else {
						unlink(c(file.path(path, 'size.txt'), 'NUL'))
						stop('Select another annotator set.')
					}
				}
		},
		both = {
			f1 <- paste0(path, '/gnomad.exomes.v4.0.sites.', chr, '.vcf.bgz')
			f2 <- paste0(path, '/gnomad.genomes.v4.0.sites.', chr, '.vcf.bgz')
			if(file.exists(f1)) {
				return
				} else  {
					lnk <- paste0(' https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/exomes/gnomad.exomes.v4.0.sites.', 
						chr, '.vcf.bgz')
					system(paste0("curl -sI", lnk, " | grep -i Content-Length | cut -d ' ' -f 2 | uniq > ", path, "/size.txt"))
					message('You are about to download ~', round(as.integer(readLines(file.path(path, 'size.txt')))*10^-9, 2), 'GiB of data.') 
					if(isTRUE(askYesNo('Do you want to proceed?'))){
						system(paste0('wget -q --show-progress -P ', path, lnk))
						lnkIndex <- paste0(lnk, '.tbi')
						system(paste0('wget -q --show-progress -P ', path, lnkIndex))
						unlink(c(file.path(path, 'size.txt'), 'NUL'))
					} else {
						unlink(c(file.path(path, 'size.txt'), 'NUL'))
						stop('Select another annotator set.')
					}
				}
			if(file.exists(f2)) {
				return
				} else  {	
					lnk <- paste0(' https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/genomes/gnomad.genomes.v4.0.sites.', 
						chr, '.vcf.bgz')	
					system(paste0("curl -sI", lnk, " | grep -i Content-Length | cut -d ' ' -f 2 | uniq > ", path, "/size.txt"))
					message('You are about to download ~', round(as.integer(readLines(file.path(path, 'size.txt')))*10^-9, 2), 'GiB of data.') 
					if(isTRUE(askYesNo('Do you want to proceed?'))){
						system(paste0('wget -q --show-progress -P ', path, lnk))
						lnkIndex <- paste0(lnk, '.tbi')
						system(paste0('wget -q --show-progress -P ', path, lnkIndex))
						unlink(c(file.path(path, 'size.txt'), 'NUL'))
					} else {	
						unlink(c(file.path(path, 'size.txt'), 'NUL'))
						stop('Select another annotator set.')
					}
				}
		}
	)
}

#' Download ClinVar files
#'
#' This function downloads ClinVar vcf file.
#'
#' @param path Where to download data
#' @return ClinVar vcf file.
#' @export
download_clinvar <- function(path){
	suppressWarnings(dir.create(path, recursive = TRUE))
	if(file.exists(file.path(path, 'clinvar.vcf.gz'))){
		return
	} else {
		lnk <- ' https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz'
		system(paste0("curl -sI", lnk, " | grep -i Content-Length | cut -d ' ' -f 2 | head -1 > ", path, "/size.txt"))
		message('You are about to download ~', round(as.integer(readLines(file.path(path, 'size.txt')))*10^-9, 2), 'GiB of data.') 
		if(isTRUE(askYesNo('Do you want to proceed?'))){
			system(paste0('wget -q --show-progress -P ', path, lnk))
			# system('mv ./data/clinvar/clinvar.vcf.gz ./data/clinvar/clinvar_`date +"%Y-%m-%d"`.vcf.gz')
			unlink(c(file.path(path, 'size.txt'), 'NUL'))
			} else {
				unlink(c(file.path(path, 'size.txt'), 'NUL'))
				stop('Select another annotator set.')
			}
	}
}

#' Download AlphaMissense tsv files
#'
#' This function downloads AlphaMissense tsv file.
#'
#' @param path Where to download data
#' @return AlphaMissense vcf files
#' @export
download_am <- function(path){
	suppressWarnings(dir.create(path, recursive = TRUE))
	if(file.exists(file.path(path, 'AlphaMissense_hg38.tsv.gz'))){
		return
	} else {
		lnk <- ' https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz'
		system(paste0("curl -sI", lnk, " | grep -i Content-Length | cut -d ' ' -f 2 | head -1 > ", path, "/size.txt"))
		message('You are about to download ~', round(as.integer(readLines(file.path(path, 'size.txt')))*10^-9, 2), 'GiB of data.') 
		if(isTRUE(askYesNo('Do you want to proceed?'))){
			system(paste0('wget -q --show-progress -P ', path, lnk))
			unlink(c(file.path(path, 'size.txt'), 'NUL'))
			} else {
				unlink(c(file.path(path, 'size.txt'), 'NUL'))
				stop('Select another annotator set.')
			}
	}
}

#' Download REVEL chromosome files.
#'
#' This function downloads REVEL vcf files per chromosome.
#'
#' @param path Where to download data
#' @param chr Chromosome name
#' @return REVEL vcf files.
#' @export
download_revel <- function(chr, path) {
	suppressWarnings(dir.create(path, recursive = TRUE))
	name <- paste0('revel-v1.3_segments_chrom_', chr, '.zip')
	lnk <- paste0(' https://zenodo.org/records/7072866/files/', name, '?download=1')
	system(paste0("curl -sI", lnk, " | grep -i Content-Length | cut -d ' ' -f 2 | head -1 > ", path, "/size.txt"))
	message('You are about to download ~', round(as.integer(readLines(file.path(path, 'size.txt')))*10^-9, 2), 'GiB of data.') 
	if(isTRUE(askYesNo('Do you want to proceed?'))){
		system(paste0('wget -q --show-progress -O ', file.path(path, name), lnk))
		unlink(c(file.path(path, 'size.txt'), 'NUL'))
		system(paste('unzip', file.path(path, name), '-d', path))
		unlink(list.files(path, pattern = 'zip', full.names = TRUE))
		unlink(c(file.path(path, 'size.txt'), 'NUL'))
		} else {
			unlink(c(file.path(path, 'size.txt'), 'NUL'))
			stop('Select another annotator set.')
		}
	return(gsub('.zip', '', file.path(path, name)))
}

#' Download LOVD3 gene variants.
#'
#' This function downloads LOVD3 gene.
#'
#' @param path Where to download data
#' @param gns Gene names
#' @return LOVD3 gene files.
#' @export
download_lovd3 <- function(gns, path){
	suppressWarnings(dir.create(path, recursive = TRUE))
	lnk <- 'https://databases.lovd.nl/shared/download/all/gene/'
	for(i in gns){
		system(paste0('wget -q --show-progress -O ', path, '/LOVD_full_download_', i, '_`date +"%Y-%m-%d"`.txt ', lnk, i))
	}
	ok <- gsub('LOVD_full_download_|_.*-.*', '', list.files(path))
	message(setdiff(gns, ok))
}