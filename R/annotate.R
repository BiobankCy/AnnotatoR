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
#' @return Annotated variants
#' @export
annotateInterVar <- function(vcf){
	p <- progressr::progressor(steps = nrow(vcf))
	out <- list()
	for(i in 1:nrow(vcf)){
		p(message = sprintf("Adding %g", i))
		tmp <- suppressMessages(httr::content(httr::GET(paste0(
			'http://wintervar.wglab.org/api_new.php?queryType=position&chr=', gsub('chr', '', vcf[i, 'CHR']), 
			'&pos=', vcf[i, 'POS'],
			'&ref=', vcf[i, 'REF'],
			'&alt=', vcf[i, 'ALT'],
			'&build=hg38'
		)), 'text'))

		if(tmp == ''){
			out[[i]] <- NA	
		} else {
		 	out[[i]] <- process_json(tmp) 
		}
	}
	return(out)
}
