#' Function to prepare drugdriveR data in appropriate fashion
#'
#' param data The input dataset, consisting of patient ID, readmission status, covariates, and medications
#' return Saves one data.table in long format
#' @examples
#' prepdata()
#'

prepdata <- function(){
########## Load in the required data (HF readmission and patient drugs)
#visit_sample  <- fread("~/Projects/HF-ml/data/second/Visit - Drugs.csv")
#hfencounters <- fread("~/Projects/HF-ml/data/second/AllHF_2010_2015.csv")

visit_sample  <- fread("~/Projects/HF-ml/data/2015-12-08/Drugs.csv")
hfencounters <- fread("~/Projects/HF-ml/data/2015-12-08/AIIHF-2010-2015.csv")


#################### Rename of columns in the datasets #####################
## Change visitdrugs variable names
#setnames(visit_sample, c("MEDICAL_RECORD_NUMBER","ENCOUNTER_VISIT_ID"), #old names
 #        c("mrn","Encounter.Number")) #new names

setnames(visit_sample, c("MASKED_MRN","MASKED_VISIT_ID"), #old names
         c("mrn","Encounter.Number")) #new names

# Change hfencounters variable names
setnames(hfencounters,c('MASKED_VISIT_ID','DISCHARGE_DATE','MASKED_MRN','DBS3_DESC','DISCHARGE_DISPOSITION','DISCHARGE_DEPARTMENT_DESC','MSDRG_CODE','MSDRG_DESCRIPTION','PRIN_DIAG_CODE','PRIN_DIAG_DESC','SEX_DESC','RACE_DESC','ZIP_CODE','ACTUAL_LENGTH_OF_STAY','UHC_INDEX_ENC_FLAG','AGE_IN_YEARS'),
         c('Encounter.Number','Discharge.Date','mrn','DBS3.Desc','Discharge.Disposition','Discharge.Department.Desc','MSDRG.Code','MSDRG.Description','Prin.Diag.Code','Prin.Diag.Desc.','sex','race','zip','Actual.Length.of.Stay','readmit','age'))


################################################################################
### Parsing drugs from input files  ############################################


#####################
##################### Extract the first part of all of the drug names,
##################### splitting at the first number
#####################

##################### Extract the first part of all of the drug names, splitting at the first number #####################
medpattern <- '^.\\w*' # Matches the first word in the string
med_results <- regexpr(pattern=medpattern, visit_sample$MEDICATION_NAME, perl=TRUE) # compute the reg. expressions
mednames <- rep(NA, length(med_results)) # create an empty vector of appropriate length to store matches
mednames <- regmatches(visit_sample$MEDICATION_NAME, med_results) # store the matched names in vector
# mednames[med_results!=1] <- regmatches(visit_sample$MEDICATION_NAME, med_results) #only needed if some names do not match pattern

visit_sample$med_name <- mednames # add vector of medical names to visit_sample dataset

# Compute number of times each MEDICATION_NAME given (in total)
visit_sample[,N:=.N, by=.(MEDICATION_NAME)]

##### Check to see output of Regular expression
setkey(visit_sample, MEDICATION_NAME)
drugnamedf <- unique(visit_sample) #extract out columns with unique drugnames from visit_sample
drugnamedf <- drugnamedf[,.(MEDICATION_NAME,med_name,N)] # Keep only MEDICATION_NAME and med_name columns
colnames(drugnamedf) <- c("EHR Medication Name", "RegEXP Medication Name","N (EHR Medication Name)")

# kable(drugnamedf)

visit_sample[,N:=NULL, by=.(MEDICATION_NAME)] # Delete the N column, which was created only for the above table
rm(medpattern, mednames, drugnamedf, med_results); # cleanup extra vector

setkey(visit_sample, mrn, med_name) # Set data.table key for faster operations
visit_sample[,nTimes:=.N, by=.(mrn,med_name)] # Compute number of times each med_name given to each MRN

#####################
##################### Extract the drug dosages from the dataset
#####################

##################### Add the dosages to datafile, using regular expressions #####################
dosepattern <- '([0-9]+\\.*[0-9]*)' # matches n.n or n
dose_results <- regexpr(pattern=dosepattern, visit_sample$MEDICATION_NAME , perl=TRUE) #compute the regexes
dosages <- rep(NA, length(dose_results)) # NA vector of appropriate length, because R just skips no-matches and returns a shorter vector
dosages[dose_results!=-1] <- regmatches(visit_sample$MEDICATION_NAME, dose_results) # input the actual match strings, leaving NAs at no-match spots
visit_sample$dose <- dosages # put our results back into the dataset
visit_sample$dose <- as.numeric(visit_sample$dose) # make the vector numeric, instead of a character vector

#####################
##################### Merge HF Readmissions and Drugs Prescribed datasets into admit.drug
#####################

hfsubset <- hfencounters[,.(Encounter.Number,Discharge.Date,mrn,Actual.Length.of.Stay,readmit,sex,race,age,zip)] # extract columns
hfsubset$Encounter.Number <- as.character(hfsubset$Encounter.Number) # Change encounter numbers to strings

# Change MRN types to be the same in both datasets
hfsubset$mrn <- as.character(hfsubset$mrn)
visit_sample$mrn <- as.character(visit_sample$mrn)

### Create new dataset admit.drug by merging hfsubset and visit sample
### This dataset is the basis for the rest of the analysis
admit.drug <- merge(visit_sample, hfsubset, by=c('mrn','Encounter.Number'))

rm(visit_sample, hfsubset, hfencounters) #cleanup old data tables

##################### Cleanup of merged dataset admit.drug #####################
admit.drug[,DOSE:=NULL] # Get rid of useless dosage instruction column
setnames(admit.drug,c('Actual.Length.of.Stay', '   '),c('LOS', 'X')) # change names to make them more wieldy
admit.drug[,X:=NULL] #get rid of useless, old index (separate from above in case we decide to keep it later)
admit.drug$readmit <- as.integer(admit.drug$readmit)
setkey(admit.drug,mrn,Encounter.Number,med_name,ORDER_DATE) #set keys on data set

### Vitally important fix: if a patient is ever readmitted, give all of their entries a "1" for readmission
admit.drug <- copy(admit.drug[,readmit:=ifelse(sum(readmit>0),1L,0L), by=.(mrn)])

################### Compute the number of medications each patient is currently taking ###################
admit.drug$med_name <- as.factor(admit.drug$med_name)

setkey(admit.drug, mrn)
admit.drug[, nDrugs:=length(unique(med_name)), by=.(mrn)]

setkey(admit.drug, mrn) # Use MRN as the key for admit.drug data table

# Clean up race-variable in dataset
admit.drug$race <- make.names(admit.drug$race)
admit.drug$race <- as.factor(admit.drug$race)
admit.drug$race <- factor(admit.drug$race, levels(admit.drug$race)[c(6,3,1,2,4,5)])

# Add sex variable as factor in dataset
admit.drug$sex <- as.factor(admit.drug$sex)

# Generate datasets with only 1 line per MRN, used for plotting
# This is because we alredy computed a group statistic (nDrugs) in the data table
# admit.drug.one <- unique(admit.drug) # get only one record per mrn (n rows == n patients)

return(admit.drug)
}