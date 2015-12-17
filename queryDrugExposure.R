#' Function that computes drug exposure enrichment with outcome
#'
#' @param dataset The input dataset, consisting of patient ID, readmission status, covariates, and medications
#' @param drugcolumnname The name of the dataset's column containing drug information
#' @param patientcolumnname The name of the dataset's column containing patient identifiers
#' @param outcomecolname The name of the dataset's column containing binary outcome
#' @param covariatecols The name of dataset columns containing covariates
#' @param useBayes Use bayesian regression from arm package (1=yes, 0=no (default) )
#' @return One xlsx file or X different CSVs with results.
#' @examples
#' queryDrugExposure(dataset=admit.drug.formula.subset,
#'                   drugcolumnname=drugcolumnname,
#'                   patientcolumnname=patientcolumnname,
#'                   outcomecolname=outcomecolname,
#'                   covariatecols=covariatecols,
#'                   useBayes=useBayes)
#'

queryDrugExposure <- function(dataset=admit.drug.formula.subset,
                              drugcolumnname=drugcolumnname,
                              patientcolumnname=patientcolumnname,
                              outcomecolname=outcomecolname,
                            covariatecols=covariatecols,
                              useBayes=useBayes
                              ){

admit.drug.exposure <- dataset
#admit.drug.exposure <- admit.drug.formula.subset

cat("Computing Drug Exposures\n");

AllDrugs <- unique(admit.drug.exposure[, get(eval(drugcolumnname))] )

cat("Constructing drug exposure variables\n")

pb <- txtProgressBar(min = 0, max = length(AllDrugs), style = 3) #initialze a progress bar

nPT <- length( unique(admit.drug.exposure[, get(eval(patientcolumnname))] ))
Exp_colnames <- paste0("Exp_", AllDrugs)
Exp_mat <- matrix(data=NA, nrow=nPT, ncol=length(Exp_colnames))
colnames(Exp_mat) <- Exp_colnames

#admit.drug.exposure <- cbind(admit.drug.exposure, Exp_mat)
#Exp_indices <- (ncol(admit.drug.exposure)-(length(Exp_colnames)-1)) : ncol(admit.drug.exposure)

#colnames(admit.drug.exposure)[min(Exp_indices):max(Exp_indices)] <- Exp_colnames
#colnames(admit.drug.exposure) <- tolower(colnames(admit.drug.exposure))

### This loop creates a 1/0 variable if the patient is receiving that drug (or not)
for(i in seq_along(AllDrugs)){
  Exp_mat[,i] <- admit.drug.exposure[, ifelse(AllDrugs[i] %in% get(eval(drugcolumnname)),1,0),by=get(eval(patientcolumnname))]$V1
  setTxtProgressBar(pb, i)
}

# system.time( # This way is much slower
# for(i in seq_along(AllDrugs)){
#  admit.drug.exposure[,eval(paste0("Exp_",AllDrugs[i])):=ifelse(AllDrugs[i] %in% get(eval(drugcolumnname)), 1, 0), by=get(eval(patientcolumnname))];
#  setTxtProgressBar(pb, i) })

close(pb)

keycols = c(patientcolumnname)
setkeyv(admit.drug.exposure,  keycols);
admit.drug.exposure.one <- unique(admit.drug.exposure); #get one row per MRN
admit.drug.exposure.one <- cbind(admit.drug.exposure.one, Exp_mat)

# Get indices of all of our drug exposure columns in the dataset
# exp_columns <- grep("Exp_", colnames(admit.drug.exposure.one), value = TRUE)

########## Loop over dataset performing regressions ###########
########## This could be parallelized if required
estimate_df <- as.data.frame(matrix(ncol=4, nrow=length(Exp_indices)), data=NA)
colnames(estimate_df) <- c("Estimate", "Std..Error", "z.value", "Pr...z..")

cat("Performing regressions\n")

pb <- txtProgressBar(min = 0, max = length(Exp_indices), style = 3)
for(i in seq_along(Exp_indices)){

#for(i in seq_along(exp_columns[1:10])){

  exposure <- parse(text=Exp_colnames[i])
  #exposure <- 'Exp_ONDANSETRON'
  regformula <- as.formula( paste(outcomecolname, "~", exposure, "+", paste(covariatecols, collapse=" + ")))

  # Do a regression, pasting in the drug exposure column as a covariate
  if(useBayes==1){
    out <- bayesglm(regformula,
                    data=admit.drug.exposure.one, family=binomial(link="logit"))
  }else{
    out <- glm(regformula,
               data=admit.drug.exposure.one, family=binomial(link="logit"))
  }

  # Append the results to our data frame
  estimate_df[i,] <- coef(summary(out))[2,];
  setTxtProgressBar(pb, i)
}
close(pb)

### Make row names and variable names
row.names(estimate_df) <- Exp_colnames
colnames(estimate_df) <- c('Estimate', 'Std.Error', "Z_value","P_value")

### Order columns by P-value
estimate_df <- estimate_df[with(estimate_df, order(P_value)), ]

### FDR significant P-values
FDR.pval <- p.adjust(estimate_df$P_value, method="fdr")
estimate_df <- cbind(estimate_df, FDR.pval)

### Rework dataset into data.table format
estimate_df <- data.table(cbind(rownames(estimate_df), estimate_df))
setnames(estimate_df, "rownames(estimate_df)", "Drug.Exposure")

### Create indicator variable for FDR Significance
estimate_df[,FDR.Sig:=ifelse(FDR.pval<0.05,1,0), by=.(Drug.Exposure)]

### Make the drug name column look better by removing "Exp_" before the drug name
estimate_df$Drug.Exposure <- gsub("^.*?_", "", as.character(estimate_df$Drug.Exposure))

### Add OR and CI to table (exponentiated coefficients)
estimate_df[, OR:=exp(Estimate)]

# Write this table to a CSV
return(estimate_df)

}

