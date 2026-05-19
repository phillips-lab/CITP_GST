# Set up ------------------------------------------------------------------

library(readxl); library(coxme); library(multcomp); 
library(lme4); library(lubridate); 

activedir <- file.choose()
setwd(dirname(activedir))

comps <- read_excel("dose_responses.xlsx")
date <- format(today(), "%m%d%y")

comps$Lab <- as.factor(comps$Lab)
comps$StartDate <- as.factor(comps$StartDate)
comps$Tech <- as.factor(comps$Tech)
comps$Compound <- as.factor(comps$Compound)
comps$Concentration <- as.factor(comps$Concentration)
comps$Strain <- as.factor(comps$Strain)
comps$death_id <- as.factor(comps$death_id)
comps$observation_id <- as.factor(comps$observation_id)
comps$plate_id <- as.factor(comps$plate_id)
comps$experiment_id <- as.factor(comps$experiment_id)
comps$Set <- as.factor(comps$Set)
comps$Comp_Conc <- droplevels(with(comps, interaction(Compound,Concentration, sep="x")))

# sets are by compound
sets <- c("setA", "setB", "setC", "setD")

strain <- c("AF16", "JU1630", "N2_PD1073")
header <- c("af16", "ju1630", "n2_pd1073")

setA <- droplevels(comps[ which(comps$Set=="GST"), ])
setB <- droplevels(comps[ which(comps$Set=="aloin"), ])
setC <- droplevels(comps[ which(comps$Set=="GLA"), ])
setD <- droplevels(comps[ which(comps$Set=="lactulose"), ])


setA_cc <- levels(setA$Comp_Conc)
setB_cc <- levels(setB$Comp_Conc)
setC_cc <- levels(setC$Comp_Conc)
setD_cc <- levels(setD$Comp_Conc)

# Pairwise Comparisons ----------------------------------------------------

jvector <- 1:4

for (j in jvector) {
  assign(paste(sets[j],"cc","pairs",sep="_"), c())
  assign(paste(sets[j],"cc","len",sep="_"), length(get(paste(sets[j],"cc",sep="_"))))
  
}

for (k in 1:(setA_cc_len-1)) {
  setA_cc_pairs[k] <- c(paste(setA_cc[setA_cc_len+1-k],"-",setA_cc[1], " = 0"))
}
for (k in 1:(setB_cc_len-1)) {
  setB_cc_pairs[k] <- c(paste(setB_cc[setB_cc_len+1-k],"-",setB_cc[1], " = 0"))
}
for (k in 1:(setC_cc_len-1)) {
  setC_cc_pairs[k] <- c(paste(setC_cc[setC_cc_len+1-k],"-",setC_cc[1], " = 0"))
}
for (k in 1:(setD_cc_len-1)) {
  setD_cc_pairs[k] <- c(paste(setD_cc[setD_cc_len+1-k],"-",setD_cc[1], " = 0"))
}

setA_cc_pairs
setB_cc_pairs
setC_cc_pairs
setD_cc_pairs


# Subset by strain ------------------------------------------------------------------

ivector <- 1:3 

for (j in jvector) {
  for (i in ivector) {
    assign(paste(sets[j],header[i],sep="_"),droplevels(get(sets[j])[ which(get(sets[j])$Strain==strain[i]), ]))
  }
}


# Cox ---------------------------------------------------------------------

for (j in jvector) {
  for (i in ivector) {
    assign(paste(sets[j],header[i],"cox",sep="_"),coxme(Surv(DeathAge,Dead) ~ Comp_Conc + (1|StartDate/Tech/plate_id), data=get(paste(sets[j],header[i],sep="_"))))
    assign(paste(sets[j],header[i],"cox","pairs",sep="_"),summary(glht(get(paste(sets[j],header[i],"cox",sep="_")), linfct = mcp(Comp_Conc = get(paste(sets[j],"cc","pairs",sep="_"))))))
  }
}


# GLM ---------------------------------------------------------------------

for (j in jvector) {
  for (i in ivector) {
    assign(paste(sets[j],header[i],"lm",sep="_"),lmer(DeathAge ~ Comp_Conc + (1|StartDate/Tech/plate_id), data=get(paste(sets[j],header[i],sep="_"))))
    assign(paste(sets[j],header[i],"lm","prof",sep="_"),profile(get(paste(sets[j],header[i],"lm",sep="_"))))
    assign(paste(sets[j],header[i],"lm","prof","CI",sep="_"),confint(get(paste(sets[j],header[i],"lm","prof",sep="_"))))
    assign(paste(sets[j],header[i],"lm","pairs",sep="_"),summary(glht(get(paste(sets[j],header[i],"lm",sep="_")), linfct = mcp(Comp_Conc = get(paste(sets[j],"cc","pairs",sep="_"))))))
  }
}


# Print to txt file -------------------------------------------------------

for (j in jvector) {
  for (i in ivector) {
    sink(paste(sets[j],"_",header[i],"_newcomps_",date,".txt",sep=""), append = FALSE)
    writeLines(paste(sets[j]," for Strain ",strain[i],sep=""))  
    writeLines("\nCox ####################\n")
    print(get(paste(sets[j],header[i],"cox",sep="_")))
    writeLines("\nCox Planned Comparisons\n")
    print(get(paste(sets[j],header[i],"cox","pairs",sep="_")))
    writeLines("\nLinear Model ####################\n")
    print(summary(get(paste(sets[j],header[i],"lm",sep="_"))))
    print(get(paste(sets[j],header[i],"lm","prof","CI",sep="_")))
    writeLines("\nLM Planned Comparisons\n")
    print(get(paste(sets[j],header[i],"lm","pairs",sep="_")))
    sink()
  }
} 


# Notes -------------------------------------------------------------------

# Quantiles and summary --------------------------------------------------------------

summary_list <- list()

for (j in seq_along(jvector)) {
  for (i in seq_along(ivector)) {

    dataset_name <- paste(sets[j], header[i], sep = "_")
    dataset <- get(dataset_name) 
    
    survfit_obj <- survfit(Surv(DeathAge, Dead) ~ Comp_Conc, data = dataset)
    assign(paste(sets[j], header[i], "survfit", sep = "_"), survfit_obj)
    
    # quantiles
    quantiles <- summary(survfit_obj)$table
    
    # calculate 90th quantile for each Comp_Conc
    quantile_90_list <- sapply(levels(dataset$Comp_Conc), function(comp_conc) {
      subset_data <- dataset[dataset$Comp_Conc == comp_conc, ]
      if (nrow(subset_data) > 0) {
        survfit_obj_comp <- survfit(Surv(DeathAge, Dead) ~ 1, data = subset_data)
        tryCatch(
          quantile(survfit_obj_comp, probs = 0.9)$quantile[1],
          error = function(e) NA
        )
      } else {
        NA
      }
    })
    
    # summary 
    summary_df <- data.frame(
      Set = sets[j],
      Strain = header[i],
      Comp_Conc = rownames(quantiles),
      n_dead = quantiles[,"events"],
      n_censored = quantiles[,"records"] - quantiles[,"events"],
      n_total = quantiles[,"records"],
      median = quantiles[,"median"],
      LCL_95 = quantiles[,"0.95LCL"],
      UCL_95 = quantiles[,"0.95UCL"],
      rmean = quantiles[,"rmean"],
      se_rmean = quantiles[,"se(rmean)"],
      quantile_90 = quantile_90_list
    )
    
    summary_list[[paste(j, i, sep = "_")]] <- summary_df
  }
}

summary_results <- do.call(rbind, summary_list)

# write the data frame to a CSV file
write.csv(summary_results, "summary_results.csv", row.names = FALSE)



