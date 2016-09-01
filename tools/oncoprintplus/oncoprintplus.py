import csv

def fetch_patient_covariates( covariate_file ):
	covariate_fh = open(covariate_file,'rb')
	covariates = csv.reader(covariate_fh, delimiter='\t')
	header = covariates.next()
	options = []
	for item in header:
		if item.endswith("_col"):
			continue
		if item == "patient":
			continue
		options.append((item,item,False))
	return options