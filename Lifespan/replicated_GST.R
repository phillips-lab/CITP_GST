# Set Up ------------------------------------------------------------------

# load packages
library(readxl); library(coxme); library(multcomp); 
library(lme4); library(lubridate); library(purrr)

# set active directory
activedir <- file.choose()
setwd(dirname(activedir))

# load data
comps <- read_excel("P3_GST.xlsx")
date <- format(today(), "%m%d%y")

# backup dataset
origcomps <- comps

strain <- c("AF16", "ED3092", "HK104", "JU775", "MY16", "N2_PD1073")
header <- c("af16", "ed3092", "hk104","ju775", "my16", "n2pd1073")

comps$Strain <- as.factor(comps$Strain)
comps$Lab <- as.factor(comps$Lab)
comps$Tech <- as.factor(comps$Tech)
comps$StartDate <- as.factor(comps$StartDate)
comps$Compound <- as.factor(comps$Compound)
comps$Concentration <- as.factor(comps$Concentration)
comps$Rep <- as.factor(comps$Rep)
comps$plate_id <- as.factor(comps$plate_id)
comps$Comp_Conc <- droplevels(with(comps, interaction(Compound,Concentration, sep="x")))

levels(comps$Strain)
levels(comps$Comp_Conc)

comps_cc_pairs <- c("gold_sodium_thiomalatex100 - CTRL_H2Ox0 = 0")

ivector <- c(1:6)

# subset data
for (i in ivector) {
  assign(paste(header[i],sep=""),droplevels(comps[ which(comps$Strain==strain[i]), ]))
}


# Cox ---------------------------------------------------------------------

for (i in ivector) {
  assign(paste(header[i],"cox",sep="_"),coxme(Surv(DeathAge,Dead) ~ Comp_Conc + (1|Lab/StartDate/Tech/plate_id), data=get(paste(header[i],sep="_"))))
  assign(paste(header[i],"cox","pairs",sep="_"),summary(glht(get(paste(header[i],"cox",sep="_")), linfct = mcp(Comp_Conc = get(paste("comps","cc","pairs",sep="_"))))))
}


# GLM ---------------------------------------------------------------------

for (i in ivector) {
  assign(paste(header[i],"lm",sep="_"),lmer(DeathAge ~ Comp_Conc + (1|Lab/StartDate/Tech/plate_id), data=get(paste(header[i],sep="_"))))
  assign(paste(header[i],"lm","prof",sep="_"),profile(get(paste(header[i],"lm",sep="_"))))
  assign(paste(header[i],"lm","prof","CI",sep="_"),confint(get(paste(header[i],"lm","prof",sep="_"))))
  assign(paste(header[i],"lm","pairs",sep="_"),summary(glht(get(paste(header[i],"lm",sep="_")), linfct = mcp(Comp_Conc = get(paste("comps","cc","pairs",sep="_"))))))
}


# text output -------------------------------------------------------------

ivector <- 1:6

for (i in ivector) {
  sink(paste(header[i],date,".txt",sep="_"), append = FALSE)
  writeLines(strain[i])  
  writeLines("\nCox\n")
  print(get(paste(header[i],"cox",sep="_")))
  print(get(paste(header[i],"cox","pairs",sep="_")))
  writeLines("\nLinear Model\n")
  print(summary(get(paste(header[i],"lm",sep="_"))))
  print(get(paste(header[i],"lm","prof","CI",sep="_")))
  print(get(paste(header[i],"lm","pairs",sep="_")))
  sink()
} 

# multiple comparisons, data table
pairsoutput <- function(modelmultcomp, model, source) {
  est   <- as.data.frame(modelmultcomp[["test"]][["coefficients"]])
  stderr<- as.data.frame(modelmultcomp[["test"]][["sigma"]])
  zval  <- as.data.frame(modelmultcomp[["test"]][["tstat"]])
  pval  <- as.data.frame(modelmultcomp[["test"]][["pvalues"]])
  
  output <- cbind(est, stderr, zval, pval, model, source)
  output$row <- row.names(output)
  output <- setNames(output, c("est","stderr","zval","pval","model","strain","conc_ctrl"))
  return(output)
}

# loop over strains automatically
strainslist<- list()

for (i in seq_along(strain)) {
  sname <- header[i]  
  cox_obj <- get(paste0(sname, "_cox_pairs"))
  lm_obj <- get(paste0(sname, "_lm_pairs"))
  df <- rbind(pairsoutput(cox_obj, "cox", sname), pairsoutput(lm_obj, "lm", sname))
  strainslist[[sname]] <- df
}

allpairs <- do.call(rbind, strainslist)

# save as excel file
write_xlsx(allpairs, paste("mult_comp_",date,".xlsx",sep=""))



# Linear Model ------------------------------------------------------------

# For overall sources of variance

lm_all <-lmer(DeathAge ~ Compound + (1|Lab/StartDate/Tech/plate_id) + (1|Species/Strain) + (1|Lab:Strain) + (1|Lab:Species) + (1|Compound:Species) + (1|Compound:Strain) + (1|Compound:Lab), data=comps)
summary(lm_all)

sink(paste("gst_variance_",date,".txt",sep=""), append = FALSE)
writeLines("\nSources of Variance - Overall\n")
summary(lm_all)
sink()


# functions for outputs  ----------------------------------------------------

# random effects variance 
lmoutput <- function(origlmmodel,model,source){
  output <- cbind(as.data.frame(summary(origlmmodel)[["varcor"]]),model,source)
  output <- setNames(output,c("group","intercept","na","variance","stddev","model","strain"))
  return(output)
}

# random effects standard deviation confidence intervals
lmcioutput <- function(origlmci,model,source){
  output <- cbind(as.data.frame(origlmci),model,source)
  output <- setNames(output,c("2.5%","97.5%","model","strain"))
  output$row <- row.names(output)
  return(output)
}

# random effects variance
coxoutput <- function(origcoxmodel,model,source){
  output <- cbind(as.data.frame(origcoxmodel[["vcoef"]]),model,source)
  output <- setNames(output,c("Lab. Date.Tech.Plate","Lab.Date.Tech","Lab.Date","Lab", "model","strain")) 
  return(output)
}


# output by function ------------------------------------------------------

# define the model object prefixes 
cox_models <- mget(paste0(header, "_cox"))
lm_models  <- mget(paste0(header, "_lm"))
lm_ci_models <- mget(paste0(header, "_lm_prof_CI"))

## variance of the random effects, cox model only
allcoxvar <- map2_dfr(cox_models, header, ~ coxoutput(.x, "cox", .y)) %>% t()
write.csv(allcoxvar, "VarianceRandomEffectsCoxPH.csv")

## variance of the random effects, linear model only-
alllmvar <- map2_dfr(lm_models, header, ~ lmoutput(.x, "lm", .y))
write.csv(alllmvar, "VarianceRandomEffectsLM.csv")

## confidence intervals of the variance of the random effects, linear model only
# only care about rows sig01 through sigma
alllmvarci <- map2_dfr(lm_ci_models, header, ~ lmcioutput(.x, "lm_CI", .y))
write.csv(alllmvarci, "VarianceRandomEffectsLMCI.csv")

