#' Review status on aggregate records (VCV and RCV) - ClinVar
revstatus_map <- data.frame(
	stars = c('four', 'three', 'two', rep('one', 2), rep('none', 3)),
	REVSTAT = c('practice_guideline', 'reviewed_by_expert_panel',
		'criteria_provided,_multiple_submitters,_no_conflicts',
		'criteria_provided,_conflicting_classifications',
		'criteria_provided,_single_submitter',
		'no_assertion_criteria_provided',
		'no_classification_provided',
		'no_classification_for_the_individual_variant')
)