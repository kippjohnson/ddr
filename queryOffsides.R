#' Function that computes offsides enrichment with outcome
#'
#' @param dataset The input dataset, consisting of patient ID, readmission status, covariates, and medications
#' @param drugcolumnname The name of the dataset's column containing drug information
#' @param patientcolumnname The name of the dataset's column containing patient identifiers
#' @param outcomecolname The name of the dataset's column containing binary outcome
#' @param covariatecols The name of dataset columns containing covariates
#' @param offsides The offides database, read into memory
#' @param useBayes Use bayesian regression from arm package (1=yes, 0=no (default) )
#' @examples
#' queryOffsides(dataset=admit.drug.formula.subset,
#'            drugcolumnname=drugcolumnname,
#'            patientcolumnname=patientcolumnname,
#'            outcomecolname=outcomecolname,
#'            covariatecols=covariatecols,
#'            offsides=offsides,
#'            useBayes=useBayes)
#'

queryOffsides <- function(dataset=admit.drug.formula.subset,
                       drugcolumnname=drugcolumnname,
                       patientcolumnname=patientcolumnname,
                       outcomecolname=outcomecolname,
                       covariatecols=covariatecols,
                       offsides=offsides,
                       useBayes=useBayes
                       ){

# admit.drug.exposure <- dataset
admit.drug.exposure <- admit.drug.formula.subset

cat("Querying Offsides for Enrichments...\n");

  # Change med_name to lower case to match lower-case names in offsides table
  #print(names(admit.drug.exposure))
  admit.drug.exposure[, med_name_lc:=NULL]
  admit.drug.exposure[, med_name_lc:=tolower( get(eval(drugcolumnname )))]

  # Set keys for admit.drug and offsides before merging
  setkey(admit.drug.exposure, med_name_lc) #set key for admit.drug as the "med_name" column, which we use to match the side effect table
  setkey(offsides, drug, event) # we'll match med_name to drug

  ###### Create subsets of the two datasets to merge together in the next step ######
  admit.drug.subset <- admit.drug.exposure[,.(mrn, med_name_lc, sex, race, age, readmit)] # extract relevant columns of admit.drug
  offsides.subset <- offsides[,.(drug, umls_id, event)] # extract relevant columns of 1-way side effects

  ############ Merge of relevant columns of admit.drug and offsides ############
  admit.drug.se <- merge(admit.drug.subset, offsides.subset, by.x='med_name_lc', by.y='drug', allow.cartesian=TRUE) #

  ############ Get rid of duplicate entries from the "cross product/cartesian" merge earlier ############
  admit.drug.se[,med_name_lc:=NULL];

  keycols = c(patientcolumnname, "event", "umls_id", covariatecols, outcomecolname)
  setkeyv(admit.drug.se,  keycols);
  admit.drug.se <- unique(admit.drug.se)


  rm(admit.drug.subset, offsides.subset);

  AllEvents <- unique(admit.drug.se$event);

  cat("Constructing event variables\n")
  pb <- txtProgressBar(min = 0, max = length(AllEvents), style = 3)

  nPT <- length( unique(admit.drug.se[, get(eval(patientcolumnname))] ))
  nDrug_colnames <- paste0("nDrugs_", AllEvents)
  nDrug_mat <- matrix(data=NA, nrow=nPT, ncol=length(nDrug_colnames))
  colnames(nDrug_mat) <- nDrug_colnames

  ## For every event, get number of drugs a patient is on with that event as a primary side effect
  for(i in seq_along(AllEvents)){
    nDrug_mat[,i] <- admit.drug.se[, ifelse(AllEvents[i] %in% event,1 ,0), by=get(eval(patientcolumnname))]$V1
    setTxtProgressBar(pb, i)
  }

#
#   ## This loop is very slow
#   ## For every event, get number of drugs a patient is on with that event as a primary side effect
#   for(i in seq_along(AllEvents)){
#     admit.drug.se[, eval(paste0("nDrugs_",AllEvents[i])):=ifelse(AllEvents[i] %in% event, 1, 0), by=eval(patientcolumnname)]
#     setTxtProgressBar(pb, i)
#   }

  close(pb)

  keycols = c(patientcolumnname)
  setkeyv(admit.drug.se,  keycols);
  admit.drug.se.one <- unique(admit.drug.se); #one line per MRN (patientcolumnname)
  admit.drug.se.one <- cbind(admit.drug.se.one, nDrug_mat)

#   # Convert regression covariates to factors (should already have been done)
#   admit.drug.se.one$race <- make.names(admit.drug.se.one$race)
#   admit.drug.se.one$race <- as.factor(admit.drug.se.one$race)
#   admit.drug.se.one$race <- factor(admit.drug.se.one$race, levels(admit.drug.se.one$race)[c(6,3,1,2,4,5)])
#   admit.drug.se.one$sex <- as.factor(admit.drug.se.one$sex)

  colnames(admit.drug.se.one) <- make.names(colnames(admit.drug.se.one))
  nDrug_columns <- grep("nDrugs_", colnames(admit.drug.se.one), value = TRUE)

  estimate_df2 <- as.data.frame(matrix(ncol=4, nrow=length(nDrug_columns)), data=NA)
  colnames(estimate_df2) <- c("Estimate", "Std..Error", "z.value", "Pr...z..")

  cat("Performing regressions\n")

  pb <- txtProgressBar(min = 0, max = length(nDrug_columns), style = 3)
  for(i in seq_along(nDrug_columns)){

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
    estimate_df2[i,] <- data.frame(coef(summary(out)))[2,]
    setTxtProgressBar(pb, i)
  }
  close(pb)

  ### Make row names variable names
  row.names(estimate_df2) <- nDrug_columns
  colnames(estimate_df2) <- c('Estimate', 'Std.Error', "Z_value","P_value")

  ### Order columns by P-value
  estimate_df2 <- estimate_df2[with(estimate_df2, order(P_value)), ]

  ### FDR significant P-values
  FDR.pval <- p.adjust(estimate_df2$P_value, method="fdr")
  estimate_df2 <- cbind(estimate_df2, FDR.pval)

  ### Rework dataset into data.table formate
  estimate_df2 <- data.table(cbind(rownames(estimate_df2), estimate_df2))
  setnames(estimate_df2, "rownames(estimate_df2)", "Effect.Exposure")

  ### Create indicator variable for FDR Significance
  estimate_df2[,FDR.Sig:=ifelse(FDR.pval<0.05,1,0), by=.(Effect.Exposure)]

  ### Make the drug name column look better by removing "nDrugs_" before the drug name
  estimate_df2$Effect.Exposure <- gsub("^.*?_", "", as.character(estimate_df2$Effect.Exposure))

  ### Add OR to table (exponentiated coefficients)
  estimate_df2[, OR:=exp(Estimate), by=.(Effect.Exposure)]

  return(estimate_df2)

  ### Table of primary side effect exposure odds ratios
  #setkey(estimate_df2, P_value) ;
  #estimate_df2[order(P_value)]  ;
  #kable(estimate_df2[P_value<0.05]) ;

}
