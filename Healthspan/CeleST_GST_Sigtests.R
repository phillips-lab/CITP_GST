# Set Up ------------------------------------------------------------------

library(gdata) # load data
library(car)
library(lme4)
library(multcomp)

activedir <- file.choose()
setwd(dirname(activedir))

comps <- read.csv("GST_CeleST_Output.csv")
origcomps <- comps

comps$Replicate <- as.factor(comps$Replicate)
comps$Age <- as.factor(comps$Age)
comps$Worm <- as.factor(comps$Worm)
comps$Comp_Age <- droplevels(with(comps, interaction(Compound,Age, sep="x")))

strain <- c("AF16","ED3092","HK104","JU775","MY16","N2")
header <- c("af16","ed3092","hk104","ju775","my16","n2pd1073")
compound <- c("DMSO","GST")
comp <- c("dmso","gst")


levels(comps$Age)
levels(comps$Comp_Age)


# Subset + Adjust DMSO -----------------------------------------------------

ivector <- 1:6
jvector <- 1:2

for (i in ivector) {
  assign(paste("comps",header[i],sep="_"),droplevels(comps[ which(comps$Strain==strain[i]), ]))
}

# LMER and ANOVA ----------------------------------------------------------

# lmer model, anova, and extra lmer for planned comparisions only (need combo instead of interaction)
# Anova function, car package. Not the same as the anova or aov functions. Must be type 3.

# af16
af16_lm <- lmer(AdjSwimmingScore ~ Compound * Age + (1|Tech/EggLay/Video), data=comps_af16, contrasts=list(Age=contr.sum, Compound=contr.sum))
af16_lm_Anova <- Anova(af16_lm, type=3)
af16_lmX <- lmer(AdjSwimmingScore ~ Comp_Age + (1|Tech/EggLay/Video), data=comps_af16, contrasts=list(Comp_Age=contr.sum))
af16_lmX_pairs <- glht(af16_lmX, linfct =  mcp(Comp_Age = c("gold_sodium_thiomalatex8 - CTRL_H2Ox8 = 0","gold_sodium_thiomalatex16 - CTRL_H2Ox16 = 0")))

# ed3092
ed3092_lm <- lmer(AdjSwimmingScore ~ Compound * Age + (1|Tech/EggLay/Video), data=comps_ed3092, contrasts=list(Age=contr.sum, Compound=contr.sum))
ed3092_lm_Anova <- Anova(ed3092_lm, type=3)
ed3092_lmX <- lmer(AdjSwimmingScore ~ Comp_Age + (1|Tech/EggLay/Video), data=comps_ed3092, contrasts=list(Comp_Age=contr.sum))
ed3092_lmX_pairs <- glht(ed3092_lmX, linfct =  mcp(Comp_Age = c("gold_sodium_thiomalatex8 - CTRL_H2Ox8 = 0","gold_sodium_thiomalatex16 - CTRL_H2Ox16 = 0")))

# hk104
hk104_lm <- lmer(AdjSwimmingScore ~ Compound * Age + (1|Tech/EggLay/Video), data=comps_hk104, contrasts=list(Age=contr.sum, Compound=contr.sum))
hk104_lm_Anova <- Anova(hk104_lm, type=3)
hk104_lmX <- lmer(AdjSwimmingScore ~ Comp_Age + (1|Tech/EggLay/Video), data=comps_hk104, contrasts=list(Comp_Age=contr.sum))
hk104_lmX_pairs <- glht(hk104_lmX, linfct =  mcp(Comp_Age = c("gold_sodium_thiomalatex8 - CTRL_H2Ox8 = 0","gold_sodium_thiomalatex16 - CTRL_H2Ox16 = 0")))


# ju775
ju775_lm <- lmer(AdjSwimmingScore ~ Compound * Age + (1|Tech/EggLay/Video), data=comps_ju775, contrasts=list(Age=contr.sum, Compound=contr.sum))
ju775_lm_Anova <- Anova(ju775_lm, type=3)
ju775_lmX <- lmer(AdjSwimmingScore ~ Comp_Age + (1|Tech/EggLay/Video), data=comps_ju775, contrasts=list(Comp_Age=contr.sum))
ju775_lmX_pairs <- glht(ju775_lmX, linfct =  mcp(Comp_Age = c("gold_sodium_thiomalatex6 - CTRL_H2Ox6 = 0","gold_sodium_thiomalatex12 - CTRL_H2Ox12 = 0")))

# my16
my16_lm <- lmer(AdjSwimmingScore ~ Compound * Age + (1|Tech/EggLay/Video), data=comps_my16, contrasts=list(Age=contr.sum, Compound=contr.sum))
my16_lm_Anova <- Anova(my16_lm, type=3)
my16_lmX <- lmer(AdjSwimmingScore ~ Comp_Age + (1|Tech/EggLay/Video), data=comps_my16, contrasts=list(Comp_Age=contr.sum))
my16_lmX_pairs <- glht(my16_lmX, linfct =  mcp(Comp_Age = c("gold_sodium_thiomalatex6 - CTRL_H2Ox6 = 0","gold_sodium_thiomalatex12 - CTRL_H2Ox12 = 0")))

# n2pd1073
n2pd1073_lm <- lmer(AdjSwimmingScore ~ Compound * Age + (1|Tech/EggLay/Video), data=comps_n2pd1073, contrasts=list(Age=contr.sum, Compound=contr.sum))
n2pd1073_lm_Anova <- Anova(n2pd1073_lm, type=3)
n2pd1073_lmX <- lmer(AdjSwimmingScore ~ Comp_Age + (1|Tech/EggLay/Video), data=comps_n2pd1073, contrasts=list(Comp_Age=contr.sum))
n2pd1073_lmX_pairs <- glht(n2pd1073_lmX, linfct =  mcp(Comp_Age = c("gold_sodium_thiomalatex6 - CTRL_H2Ox6 = 0","gold_sodium_thiomalatex12 - CTRL_H2Ox12 = 0")))


# LMM Confidence Intervals ------------------------------------------------

af16_lm_prof_CI <- confint(profile(af16_lm))
ed3092_lm_prof_CI <- confint(profile(ed3092_lm))
hk104_lm_prof_CI <- confint(profile(hk104_lm))
ju775_lm_prof_CI <- confint(profile(ju775_lm))
my16_lm_prof_CI <- confint(profile(my16_lm))
n2pd1073_lm_prof_CI <- confint(profile(n2pd1073_lm))


# Outputs -----------------------------------------------------------------

sink("gst_CeleST_AdjSwimmingScore_ANOVA.txt")
writeLines("AF16")
af16_lm_Anova
writeLines("")
writeLines("ED3092")
ed3092_lm_Anova
writeLines("")
writeLines("HK104")
hk104_lm_Anova
writeLines("")
writeLines("JU775")
ju775_lm_Anova
writeLines("")
writeLines("MY16")
my16_lm_Anova
writeLines("")
writeLines("N2PD1073")
n2pd1073_lm_Anova
sink()

sink("gst_CeleST_AdjSwimmingScore_LMM.txt")
writeLines("AF16")
summary(af16_lm)
writeLines("ED3092")
summary(ed3092_lm)
writeLines("HK104")
summary(hk104_lm)
writeLines("JU775")
summary(ju775_lm)
writeLines("MY16")
summary(my16_lm)
writeLines("N2PD1073")
summary(n2pd1073_lm)
sink()

sink("gst_CeleST_AdjSwimmingScore_LMM_pairs.txt")
writeLines("AF16")
summary(af16_lmX_pairs)
writeLines("ED3092")
summary(ed3092_lmX_pairs)
writeLines("HK104")
summary(hk104_lmX_pairs)
writeLines("JU775")
summary(ju775_lmX_pairs)
writeLines("MY16")
summary(my16_lmX_pairs)
writeLines("N2PD1073")
summary(n2pd1073_lmX_pairs)
sink()

sink("gst_CeleST_AdjSwimmingScore_LMM_CI.txt")
writeLines("AF16")
af16_lm_prof_CI
writeLines("ED3092")
ed3092_lm_prof_CI
writeLines("HK104")
hk104_lm_prof_CI
writeLines("JU775")
ju775_lm_prof_CI
writeLines("MY16")
my16_lm_prof_CI
writeLines("N2PD1073")
n2pd1073_lm_prof_CI
sink()


# for each
summary(af16_lm)
summary(af16_lmX)
summary(af16_lmX_pairs)
af16_lm_Anova



# summary -----------------------------------------------------------------

library(dplyr)

# Create summary table
summary_df <- comps %>%
  group_by(Strain, Compound, Age) %>%
  summarise(
    total = n(),
    mean_SwimmingScore = mean(SwimmingScore, na.rm = TRUE),
    se_SwimmingScore   = sd(SwimmingScore, na.rm = TRUE) / sqrt(sum(!is.na(SwimmingScore))),
    mean_AdjSwimmingScore = mean(AdjSwimmingScore, na.rm = TRUE),
    se_AdjSwimmingScore   = sd(AdjSwimmingScore, na.rm = TRUE) / sqrt(sum(!is.na(AdjSwimmingScore))),
    .groups = "drop"
  )

write.csv(summary_df, "SwimmingScore_Summary.csv", row.names = FALSE)


