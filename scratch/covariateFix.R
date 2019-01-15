# The GTEX-supplied covariates are all only subject ID, whereas I want them all to have SAMPLE IDs. I'm going to fix this
# by checking to see if the tissue type found in the covariate file name and the subject ID match in each row of tissue_table.txt. If there
# is a match, the column name in the GTEX covariate file will be renamed to the sample ID. This /should/ happen for all entries.
require("data.table")
require("R.utils")
args = commandArgs(trailingOnly=TRUE)
# args[1] is each covariate file launched as part of a batch script
GTEXcov <- fread("GTEx_Analysis_v7_eQTL_covariates/Adipose_Visceral_Omentum.v7.covariates.txt")
tisstabl <- fread("tissue_table.txt")

ind <- match(names(GTEXcov), tisstabl$submitted_subject_id)
names(GTEXcov) <- tisstabl$Sample_Name[ind]