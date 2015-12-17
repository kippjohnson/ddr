#' Top-level example script
#' Kipp Johnson
#' kipp.johnson@icahn.mssm.edu
#'

library(data.table, quietly = TRUE)
library(arm, quietly = TRUE)

### Configuration variables
#datapath='~/Projects/HF-ml/data/admitdrug.csv'

# Input file paths
datapath='~/Projects/HF-ml/data/2015-12-08/admitdrug.csv' # take two
outputdir='~/Desktop/'
siderpath="~/Projects/HF-ml/data/sider-full.csv"
offsidespath='~/Projects/HF-ml/data/3003377s-offsides.tsv'
twosidespath='~/Projects/HF-ml/data/3003377s-twosides.tsv'

# Data information
drugcolumnname='med_name'
patientcolumnname='mrn'
outcomecolname="readmit"
covariatecols=c("age","sex","race")

# Options
useBayes=0

  ### Input data preparation ###
  admit.drug <- fread(datapath, header=TRUE)
  xnames <- c(patientcolumnname, drugcolumnname, covariatecols)
  admit.drug.formula.subset <- admit.drug[, c(outcomecolname, xnames), with=FALSE]

  ### Database locations ####
  cat("Reading in Sider...\n"); sider <- fread(siderpath)
  #cat("Reading in Offsides...\n"); offsides <- fread(offsidespath, sep="\t")
  #cat("Reading in Twosides...\n"); suppressWarnings( twosides <- fread(twosidespath, sep="\t") )

  ### Storage flags ###
  drug_exposure_result <- paste0(outputdir, "drug_enrichment.csv")
  sider_side_effect <- paste0(outputdir, "sider_se_enrichment.csv")
  primary_side_effect <- paste0(outputdir, "primary_se_enrichment.csv")
  secondary_side_effect <- paste0(outputdir, "secondary_se_enrichment.csv")

  ### Calculate drug enrichments ###
  drug_exposure <- queryDrugExposure(dataset=admit.drug.formula.subset,
                                     drugcolumnname=drugcolumnname,
                                     patientcolumnname=patientcolumnname,
                                     outcomecolname=outcomecolname,
                                     covariatecols=covariatecols,
                                     useBayes=useBayes);
#
#   ### Calculate sider side effect enrichments ###
  sider_out <- querySider(dataset=admit.drug.formula.subset,
                          drugcolumnname=drugcolumnname,
                          patientcolumnname=patientcolumnname,
                          outcomecolname=outcomecolname,
                          covariatecols=covariatecols,
                          sider=sider,
                          useBayes=useBayes);

  ### Calculate offsides side effect enrichments ###
#   offsides_out <- queryOffsides(dataset=admit.drug.formula.subset,
#                                 drugcolumnname=drugcolumnname,
#                                 patientcolumnname=patientcolumnname,
#                                 outcomecolname=outcomecolname,
#                                 covariatecols=covariatecols,
#                                 offsides=offsides,
#                                 useBayes=useBayes);

  ### Calculate twosides side effect enrichments ###
#   twosides_out <- function(dataset=admit.drug.formula.subset,
#                             drugcolumnname=drugcolumnname,
#                             patientcolumnname=patientcolumnname,
#                             outcomecolname=outcomecolname,
#                             covariatecols=covariatecols,
#                             twosides=twosides,
#                             useBayes=useBayes)

  write.csv(drug_exposure, file=drug_exposure_result)
  write.csv(sider_out, file=sider_side_effect)
  write.csv(offsides_out, file=primary_side_effect)
  write.csv(twosider_out, file=secondary_side_effect)

