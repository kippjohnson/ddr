#' Function that computes twosides enrichment with outcome
#'
#' @param dataset The input dataset, consisting of patient ID, readmission status, covariates, and medications
#' @param drugcolumnname The name of the dataset's column containing drug information
#' @param patientcolumnname The name of the dataset's column containing patient identifiers
#' @param outcomecolname The name of the dataset's column containing binary outcome
#' @param covariatecols The name of dataset columns containing covariates
#' @param twosides The twosides database, read into memory
#' @param useBayes Use bayesian regression from arm package (1=yes, 0=no (default) )
#' @examples
#' queryTwosides(dataset=admit.drug.formula.subset,
#'               drugcolumnname=drugcolumnname,
#'               patientcolumnname=patientcolumnname,
#'               outcomecolname=outcomecolname,
#'               covariatecols=covariatecols,
#'               twosides=twosides,
#'               useBayes=useBayes)

queryTwosides <- function(dataset=admit.drug.formula.subset,
                          drugcolumnname=drugcolumnname,
                          patientcolumnname=patientcolumnname,
                          outcomecolname=outcomecolname,
                          covariatecols=covariatecols,
                          twosides=twosides,
                          useBayes=useBayes
){

# admit.drug.subset <- dataset
admit.drug.subset <- admit.drug.formula.subset

# Change med_name to lower case to match lower-case names in offsides table
#print(names(admit.drug.exposure))
admit.drug.subset[, med_name_lc:=NULL]
admit.drug.subset[, med_name_lc:=tolower( get(eval(drugcolumnname )))]

# #################################################################
# ##### Looking for intersection of drugnames between our HF and side effect datasets #####
# admit.drug.drugnames <- (unique(admit.drug.subset$med_name_lc)); # HF dataset
# twosides.drugnames <- unique(c(twosides$drug1, twosides$drug2)); # twosides side effect dataset
# intersect.drugnames <- intersect( admit.drug.drugnames, twosides.drugnames ) # intersection of drug names
#
# kable(data.frame('Dataset Drug Count'=c(length(admit.drug.drugnames),
#                                         length(twosides.drugnames),
#                                         length(intersect.drugnames)),
#                  row.names = c("HF Admissions","Side Effects","Intersection")))
#
# rm(admit.drug.drugnames, twosides.drugnames, intersect.drugnames);
# #################################################################

twosides.subset <- twosides[,.(drug1, drug2, event_name)]
setkey(twosides.subset, drug1, drug2);

# Create a new columns, with 1/0 if MRNs have >1 drug per MRN
admit.drug.subset[, inc:=ifelse(length(unique(med_name_lc))>1,1,0), by=get(eval(patientcolumnname))]

# Extract only the MRNs which have >1 drug
admit.drug.subset <- subset(admit.drug.subset, inc==1)

### Create columns of all drug-drug combinations for the _first_ patient
df_subset <- subset(admit.drug.subset, mrn==unique(admit.drug.subset$mrn)[1])

df <- (cbind(df_subset$mrn, df_subset$readmit, df_subset$race, df_subset$age, df_subset$sex, as.vector(combn(df_subset$med_name_lc, 2, FUN=function(x){c(paste0(x[1],x[2]),paste0(x[2],x[1]))})) ))


for(i in 2:length(unique(admit.drug.subset$mrn))){
  df_subset <- subset(admit.drug.subset, mrn==unique(admit.drug.subset$mrn)[i]) #take n MRN different subsets of the data
  df2 <- (cbind(df_subset$mrn, df_subset$readmit, df_subset$race, df_subset$age, df_subset$sex, as.vector(combn(df_subset$med_name_lc, 2, FUN=function(x){c(paste0(x[1],x[2]),paste0(x[2],x[1]))})) )) # bind the mrn, readmission, and unique drug combinations together
  df <- rbind(df, df2) #add the newly computed rows to the old data frame
}

admit.drug.combos <- as.data.table(df)
setnames(admit.drug.combos, c("V1","V2","V3","V4","V5","V6"), c("mrn","readmit", "race", "age", "sex", "drugcombo"))
setkey(admit.drug.combos, mrn, drugcombo)

#########################################################################################
#########################################################################################

twosides.subset$drugcombo <- paste0(twosides.subset$drug1, twosides.subset$drug2)
setkey(twosides.subset, drugcombo, event_name)

admit.drug.sse <- merge(twosides.subset, admit.drug.combos, by.x="drugcombo", by.y="drugcombo", allow.cartesian = TRUE)

setkey(admit.drug.sse, mrn, event_name, drugcombo, race, age, sex, readmit)

# This is optional (different question)
admit.drug.sse <- unique(admit.drug.sse)

#########################################################################################


AllEvents <- unique(admit.drug.sse$event_name);

pb <- txtProgressBar(min = 0, max = length(AllEvents), style = 3)

nPT <- length( unique(admit.drug.sse[, get(eval(patientcolumnname))] ))
nDrug_colnames <- paste0("nDrugs_", AllEvents)
nDrug_mat <- matrix(data=NA, nrow=nPT, ncol=length(nDrug_colnames))
colnames(nDrug_mat) <- nDrug_colnames

## For every event, get number of drugs a patient is on with that event as a primary side effect
for(i in seq_along(AllEvents)){
  nDrug_mat[,i] <- admit.drug.sse[, ifelse(AllEvents[i] %in% event_name,1 ,0), by=get(eval(patientcolumnname))]$V1
  setTxtProgressBar(pb, i)
}

# ## This loop is very slow
# ## For every event, get number of drugs a patient is on with that event as a primary side effect
# for(i in seq_along(AllEvents)){
#   admit.drug.sse[, eval(paste0("nCombos_",AllEvents[i])):=ifelse(AllEvents[i] %in% event_name, 1, 0), by=.(mrn)]
#   setTxtProgressBar(pb, i)
# }
close(pb)


###########################################################################################
###########################################################################################
###########################################################################################

setkey(admit.drug.sse, mrn);
admit.drug.sse.one <- unique(admit.drug.sse); #one line per MRN
admit.drug.sse.one <- cbind(admit.drug.sse.one, nDrug_mat)

# Convert regression covariates to factors
admit.drug.sse.one$race <- make.names(admit.drug.sse.one$race)
admit.drug.sse.one$race <- as.factor(admit.drug.sse.one$race)
admit.drug.sse.one$race <- factor(admit.drug.sse.one$race, levels(admit.drug.sse.one$race)[c(6,3,1,2,4,5)])

admit.drug.sse.one$readmit <- as.factor(admit.drug.sse.one$readmit)
admit.drug.sse.one$sex <- as.factor(admit.drug.sse.one$sex)
admit.drug.sse.one$age <- as.integer(admit.drug.sse.one$age)

colnames(admit.drug.sse.one) <- make.names(colnames(admit.drug.sse.one))
nCombos_columns <- grep("nDrugs_", colnames(admit.drug.sse.one), value = TRUE)

estimate_df3 <- as.data.frame(matrix(ncol=4, nrow=length(nCombos_columns)), data=NA)
colnames(estimate_df3) <- c("Estimate", "Std..Error", "z.value", "Pr...z..")

pb <- txtProgressBar(min = 0, max = length(nCombos_columns), style = 3)

for(i in seq_along(nCombos_columns)){

  exposure <- parse(text=nDrug_columns[i])
  regformula <- as.formula( paste(outcomecolname, "~", exposure, "+", paste(covariatecols, collapse=" + ")))

  # Do a regression, pasting in the drug exposure column as a covariate
  if(useBayes==1){
    out <- bayesglm(regformula,
                    data=admit.drug.se.one, family=binomial(link="logit") )
  }else{
    out <- glm(regformula,
               data=admit.drug.se.one, family=binomial(link="logit") )
  }

  # Append the results to our data frame
  estimate_df3[i,] <- data.frame(coef(summary(out)))[2,]
  setTxtProgressBar(pb, i)
}
close(pb)

###########################################################################################
###########################################################################################

### Make row names variable names
row.names(estimate_df3) <- nCombos_columns
colnames(estimate_df3) <- c('Estimate', 'Std.Error', "Z_value","P_value")

### Order columns by P-value
estimate_df3 <- estimate_df3[with(estimate_df3, order(P_value)), ]

### FDR significant P-values
FDR.pval <- p.adjust(estimate_df3$P_value, method="fdr")
estimate_df3 <- cbind(estimate_df3, FDR.pval)

### Rework dataset into data.table formate
estimate_df3 <- data.table(cbind(rownames(estimate_df3), estimate_df3))
setnames(estimate_df3, "rownames(estimate_df3)", "Sec.Effect.Exposure")

### Create indicator variable for FDR Significance
estimate_df3[,FDR.Sig:=ifelse(FDR.pval<0.05,1,0), by=.(Sec.Effect.Exposure)]

### Make the drug name column look better by removing "nDrugs_" before the drug name
estimate_df3$Sec.Effect.Exposure <- gsub("^.*?_", "", as.character(estimate_df3$Sec.Effect.Exposure))

### Add OR to table (exponentiated coefficients)
estimate_df3[, OR:=exp(Estimate), by=.(Sec.Effect.Exposure)]

return(estimate_df3)

}
