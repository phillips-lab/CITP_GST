
library(readxl) # load data
library(MASS) # lda
library(BBmisc) # normalize data

activedir <- file.choose()
setwd(dirname(activedir))

comps <- read_excel("GST_CeleST.xlsx")

strain <- c("AF16","ED3092","HK104","JU775","MY16","N2")
header <- c("af16","ed3092","hk104","ju775","my16","n2")
compound <- c("h2o","GST")
comp <- c("h2o","gst")
trait <- c("Wave","Body","Asymmetry","Stretch","Curling","Travel","Brush","Activity")


# GST  --------------------------------------------------------------------

comps <- read_excel("GST_CeleST.xlsx")
colnames(comps) <- gsub(" ", "_", colnames(comps))

comps$Egg_Lay_Date <- as.factor(comps$Egg_Lay_Date)
comps$Adult_Age <- as.factor(comps$Adult_Age)

names(comps)[names(comps) == "Egg_Lay_Date"] <- "EggLay"
names(comps)[names(comps) == "Adult_Age"] <- "Age"
names(comps)[names(comps) == "Trial"] <- "Replicate"

names(comps)[names(comps) == "Wave_Init_Rate"] <- "Wave"
names(comps)[names(comps) == "Body_Wave_Number"] <- "Body"
names(comps)[names(comps) == "Asymmetry"] <- "Asymmetry"
names(comps)[names(comps) == "Stretch"] <- "Stretch"
names(comps)[names(comps) == "Curling"] <- "Curling"
names(comps)[names(comps) == "Travel_Speed"] <- "Travel"
names(comps)[names(comps) == "Brush_Stroke"] <- "Brush"
names(comps)[names(comps) == "Activity_Index"] <- "Activity"

comps$Compound <- as.factor(comps$Compound)
comps$Worm <- as.factor(comps$Worm)
comps$Replicate <- as.factor(comps$Replicate)
comps$Strain <- as.factor(comps$Strain)
comps$Tech <- as.factor(comps$Tech)

levels(comps$Compound)
levels(comps$Tech)
levels(comps$Strain)


# Subset and Adjust for h2o (SULF) -----------------------------------------

ivector <- 1:6 #numbers of strains
jvector <- 1:2 #number of compounds

# two distinct groups: subset entire dataset for each strain (starts with "comps"), and subset just the h2o data ( starts with"h2o")
# the "h2o" headed variables are used to get the lda loadings and the h2o sd's... which are later applies to the comps_strain dataframes.

for (i in ivector) {
  assign(paste("h2o",header[i],sep="_"),droplevels(comps[ which(comps$Strain==strain[i] & comps$Compound=="CTRL_H2O"), ]))
  
  assign(paste("h2o",header[i],trait[1],"sd",sep="_"),sd(get(paste("h2o",header[i],sep="_"))$Wave, na.rm = TRUE))
  assign(paste("h2o",header[i],trait[2],"sd",sep="_"),sd(get(paste("h2o",header[i],sep="_"))$Body, na.rm = TRUE))
  assign(paste("h2o",header[i],trait[3],"sd",sep="_"),sd(get(paste("h2o",header[i],sep="_"))$Asymmetry, na.rm = TRUE))
  assign(paste("h2o",header[i],trait[4],"sd",sep="_"),sd(get(paste("h2o",header[i],sep="_"))$Stretch, na.rm = TRUE))
  assign(paste("h2o",header[i],trait[5],"sd",sep="_"),sd(get(paste("h2o",header[i],sep="_"))$Curling, na.rm = TRUE))
  assign(paste("h2o",header[i],trait[6],"sd",sep="_"),sd(get(paste("h2o",header[i],sep="_"))$Travel, na.rm = TRUE))
  assign(paste("h2o",header[i],trait[7],"sd",sep="_"),sd(get(paste("h2o",header[i],sep="_"))$Brush, na.rm = TRUE))
  assign(paste("h2o",header[i],trait[8],"sd",sep="_"),sd(get(paste("h2o",header[i],sep="_"))$Activity, na.rm = TRUE))
  
  assign(paste("comps",header[i],sep="_"),droplevels(comps[ which(comps$Strain==strain[i]), ]))
  
  assign(paste("h2o",header[i],"sd",sep="_") , c(get(paste("h2o",header[i],trait[1],"sd",sep="_")),get(paste("h2o",header[i],trait[2],"sd",sep="_")),get(paste("h2o",header[i],trait[3],"sd",sep="_")),get(paste("h2o",header[i],trait[4],"sd",sep="_")),get(paste("h2o",header[i],trait[5],"sd",sep="_")),get(paste("h2o",header[i],trait[6],"sd",sep="_")),get(paste("h2o",header[i],trait[7],"sd",sep="_")),get(paste("h2o",header[i],trait[8],"sd",sep="_"))))
  
}

# scale h2o before LDA
h2o_af16[14:21] <- t(t(h2o_af16[14:21]) / h2o_af16_sd)
h2o_ed3092[14:21] <- t(t(h2o_ed3092[14:21]) / h2o_ed3092_sd)
h2o_hk104[14:21] <- t(t(h2o_hk104[14:21]) / h2o_hk104_sd)
h2o_ju775[14:21] <- t(t(h2o_ju775[14:21]) / h2o_ju775_sd)
h2o_my16[14:21] <- t(t(h2o_my16[14:21]) / h2o_my16_sd)
h2o_n2[14:21] <- t(t(h2o_n2[14:21]) / h2o_n2_sd)


# LDA - Adjusted h2o ------------------------------------------------------

for (i in ivector) {
  assign(paste("h2o",header[i],"lda",sep="_"), lda(Age ~ Wave + Body + Asymmetry + Stretch + Curling + Travel + Brush + Activity, data=get(paste("h2o",header[i],sep="_"))))
  assign(paste("h2o",header[i],"lda","coefLD1",sep="_"),get(paste("h2o",header[i],"lda",sep="_"))$scaling[1:8,1])
  
  assign(paste("h2o",header[i],"lda","svd",sep="_"),get(paste("h2o",header[i],"lda",sep="_"))$svd)
  assign(paste("h2o",header[i],"lda","propLD1",sep="_"),(get((paste("h2o",header[i],"lda","svd",sep="_")))[1]^2) / sum((get((paste("h2o",header[i],"lda","svd",sep="_"))))^2))
  assign(paste("h2o",header[i],"means",sep="_"),get(paste("h2o",header[i],"lda",sep="_"))$means)
}

h2o_af16_lda
h2o_ed3092_lda
h2o_hk104_lda
h2o_ju775_lda
h2o_my16_lda
h2o_n2_lda


# Alter Loadings ----------------------------------------------------------

# * only if loading reversed* 

h2o_af16_lda_coefLD1 <- -h2o_af16_lda_coefLD1
h2o_ed3092_lda_coefLD1 <- -h2o_ed3092_lda_coefLD1
h2o_hk104_lda_coefLD1 <- -h2o_hk104_lda_coefLD1
h2o_ju775_lda_coefLD1 <- -h2o_ju775_lda_coefLD1
h2o_my16_lda_coefLD1 <- -h2o_my16_lda_coefLD1
h2o_n2_lda_coefLD1 <- -h2o_n2_lda_coefLD1



# Create Health Index -----------------------------------------------------

# Confirm the correct loadings
h2o_af16_lda_coefLD1
h2o_ed3092_lda_coefLD1
h2o_hk104_lda_coefLD1
h2o_ju775_lda_coefLD1
h2o_my16_lda_coefLD1
h2o_n2_lda_coefLD1

# starting with the completely unadjusted strain subsets. h2o inside of comps_strain is NOT yet adjusted.

#af16
comps_af16$SwimmingScore <- comps_af16$Wave * h2o_af16_lda_coefLD1[1] * (1/h2o_af16_sd[1]) + 
  comps_af16$Body * h2o_af16_lda_coefLD1[2] * (1/h2o_af16_sd[2]) + 
  comps_af16$Asymmetry * h2o_af16_lda_coefLD1[3] * (1/h2o_af16_sd[3]) + 
  comps_af16$Stretch * h2o_af16_lda_coefLD1[4] * (1/h2o_af16_sd[4]) + 
  comps_af16$Curling * h2o_af16_lda_coefLD1[5] * (1/h2o_af16_sd[5]) + 
  comps_af16$Travel * h2o_af16_lda_coefLD1[6] * (1/h2o_af16_sd[6]) + 
  comps_af16$Brush * h2o_af16_lda_coefLD1[7] * (1/h2o_af16_sd[7]) + 
  comps_af16$Activity * h2o_af16_lda_coefLD1[8] * (1/h2o_af16_sd[8]) 

#ed3092
comps_ed3092$SwimmingScore <- comps_ed3092$Wave * h2o_ed3092_lda_coefLD1[1] * (1/h2o_ed3092_sd[1]) + 
  comps_ed3092$Body * h2o_ed3092_lda_coefLD1[2] * (1/h2o_ed3092_sd[2]) + 
  comps_ed3092$Asymmetry * h2o_ed3092_lda_coefLD1[3] * (1/h2o_ed3092_sd[3]) + 
  comps_ed3092$Stretch * h2o_ed3092_lda_coefLD1[4] * (1/h2o_ed3092_sd[4]) + 
  comps_ed3092$Curling * h2o_ed3092_lda_coefLD1[5] * (1/h2o_ed3092_sd[5]) + 
  comps_ed3092$Travel * h2o_ed3092_lda_coefLD1[6] * (1/h2o_ed3092_sd[6]) + 
  comps_ed3092$Brush * h2o_ed3092_lda_coefLD1[7] * (1/h2o_ed3092_sd[7]) + 
  comps_ed3092$Activity * h2o_ed3092_lda_coefLD1[8] * (1/h2o_ed3092_sd[8]) 

#hk104
comps_hk104$SwimmingScore <- comps_hk104$Wave * h2o_hk104_lda_coefLD1[1] * (1/h2o_hk104_sd[1]) + 
  comps_hk104$Body * h2o_hk104_lda_coefLD1[2] * (1/h2o_hk104_sd[2]) + 
  comps_hk104$Asymmetry * h2o_hk104_lda_coefLD1[3] * (1/h2o_hk104_sd[3]) + 
  comps_hk104$Stretch * h2o_hk104_lda_coefLD1[4] * (1/h2o_hk104_sd[4]) + 
  comps_hk104$Curling * h2o_hk104_lda_coefLD1[5] * (1/h2o_hk104_sd[5]) + 
  comps_hk104$Travel * h2o_hk104_lda_coefLD1[6] * (1/h2o_hk104_sd[6]) + 
  comps_hk104$Brush * h2o_hk104_lda_coefLD1[7] * (1/h2o_hk104_sd[7]) + 
  comps_hk104$Activity * h2o_hk104_lda_coefLD1[8] * (1/h2o_hk104_sd[8]) 


#ju775
comps_ju775$SwimmingScore <- comps_ju775$Wave * h2o_ju775_lda_coefLD1[1] * (1/h2o_ju775_sd[1]) + 
  comps_ju775$Body * h2o_ju775_lda_coefLD1[2] * (1/h2o_ju775_sd[2]) + 
  comps_ju775$Asymmetry * h2o_ju775_lda_coefLD1[3] * (1/h2o_ju775_sd[3]) + 
  comps_ju775$Stretch * h2o_ju775_lda_coefLD1[4] * (1/h2o_ju775_sd[4]) + 
  comps_ju775$Curling * h2o_ju775_lda_coefLD1[5] * (1/h2o_ju775_sd[5]) + 
  comps_ju775$Travel * h2o_ju775_lda_coefLD1[6] * (1/h2o_ju775_sd[6]) + 
  comps_ju775$Brush * h2o_ju775_lda_coefLD1[7] * (1/h2o_ju775_sd[7]) + 
  comps_ju775$Activity * h2o_ju775_lda_coefLD1[8] * (1/h2o_ju775_sd[8]) 

#my16
comps_my16$SwimmingScore <- comps_my16$Wave * h2o_my16_lda_coefLD1[1] * (1/h2o_my16_sd[1]) + 
  comps_my16$Body * h2o_my16_lda_coefLD1[2] * (1/h2o_my16_sd[2]) + 
  comps_my16$Asymmetry * h2o_my16_lda_coefLD1[3] * (1/h2o_my16_sd[3]) + 
  comps_my16$Stretch * h2o_my16_lda_coefLD1[4] * (1/h2o_my16_sd[4]) + 
  comps_my16$Curling * h2o_my16_lda_coefLD1[5] * (1/h2o_my16_sd[5]) + 
  comps_my16$Travel * h2o_my16_lda_coefLD1[6] * (1/h2o_my16_sd[6]) + 
  comps_my16$Brush * h2o_my16_lda_coefLD1[7] * (1/h2o_my16_sd[7]) + 
  comps_my16$Activity * h2o_my16_lda_coefLD1[8] * (1/h2o_my16_sd[8]) 

#n2
comps_n2$SwimmingScore <- comps_n2$Wave * h2o_n2_lda_coefLD1[1] * (1/h2o_n2_sd[1]) + 
  comps_n2$Body * h2o_n2_lda_coefLD1[2] * (1/h2o_n2_sd[2]) + 
  comps_n2$Asymmetry * h2o_n2_lda_coefLD1[3] * (1/h2o_n2_sd[3]) + 
  comps_n2$Stretch * h2o_n2_lda_coefLD1[4] * (1/h2o_n2_sd[4]) + 
  comps_n2$Curling * h2o_n2_lda_coefLD1[5] * (1/h2o_n2_sd[5]) + 
  comps_n2$Travel * h2o_n2_lda_coefLD1[6] * (1/h2o_n2_sd[6]) + 
  comps_n2$Brush * h2o_n2_lda_coefLD1[7] * (1/h2o_n2_sd[7]) + 
  comps_n2$Activity * h2o_n2_lda_coefLD1[8] * (1/h2o_n2_sd[8]) 

comps_af16$SwimmingScoreMean <- mean(comps_af16$SwimmingScore,na.rm = TRUE)
comps_ed3092$SwimmingScoreMean <- mean(comps_ed3092$SwimmingScore,na.rm = TRUE)
comps_hk104$SwimmingScoreMean <- mean(comps_hk104$SwimmingScore,na.rm = TRUE)
comps_ju775$SwimmingScoreMean <- mean(comps_ju775$SwimmingScore,na.rm = TRUE)
comps_my16$SwimmingScoreMean <- mean(comps_my16$SwimmingScore,na.rm = TRUE)
comps_n2$SwimmingScoreMean <- mean(comps_n2$SwimmingScore,na.rm = TRUE)


Final_GST_Output <- rbind(comps_af16,comps_ed3092,comps_hk104,comps_ju775,comps_my16,comps_n2)

Final_GST_Output$AdjSwimmingScore <- Final_GST_Output$SwimmingScore - Final_GST_Output$SwimmingScoreMean


# Merge and Export --------------------------------------------------------

write.csv(x = Final_GST_Output, file = "GST_CeleST_Output.csv", row.names = FALSE)

# LDA Output --------------------------------------------------------------

ivector <- 1:6

sink(paste("GST_CeleST_LDA_A.txt",sep=""), append = FALSE)

for (i in ivector) {
  writeLines(paste(strain[i],sep=" "))
  writeLines("SD by Trait")
  print(get(paste("h2o",header[i],"sd",sep="_")))
  writeLines("LDA Loadings")
  print(get(paste("h2o",header[i],"lda","coefLD1",sep="_")))
  writeLines("LD1 Proportion of Trace")
  print(get(paste("h2o",header[i],"lda","propLD1",sep="_")))
  writeLines("")
}


sink()


sink(paste("GST_CeleST_LDA_B.txt",sep=""), append = FALSE)

writeLines("SD per Strain")
for (i in ivector) {
  writeLines(paste(strain[i],sep=" "))
  print(get(paste("h2o",header[i],"sd",sep="_")))
}

writeLines("\nLDA Loadings per Strain")
for (i in ivector) {
  writeLines(paste(strain[i],sep=" "))
  print(get(paste("h2o",header[i],"lda","coefLD1",sep="_")))
}

writeLines("\nLD1 Proportion of Trace per Strain")
for (i in ivector) {
  writeLines(paste(strain[i],sep=" "))
  print(get(paste("h2o",header[i],"lda","propLD1",sep="_")))
}

sink()




# MANOVA on Same LDA Set --------------------------------------------------

#af16
h2o_af16_Wave <- h2o_af16$Wave
h2o_af16_Brush <- h2o_af16$Brush
h2o_af16_Activity <- h2o_af16$Activity
h2o_af16_Travel <- h2o_af16$Travel
h2o_af16_Curling <- h2o_af16$Curling
h2o_af16_Asymmetry <- h2o_af16$Asymmetry
h2o_af16_Stretch <- h2o_af16$Stretch
h2o_af16_Body <- h2o_af16$Body

af16_manova <- manova(cbind(h2o_af16_Wave,h2o_af16_Brush,h2o_af16_Activity,h2o_af16_Travel,h2o_af16_Curling,h2o_af16_Asymmetry,h2o_af16_Stretch,h2o_af16_Body) ~ Age, data = h2o_af16)
#summary(af16_manova)

# ed3092
h2o_ed3092_Wave <- h2o_ed3092$Wave
h2o_ed3092_Brush <- h2o_ed3092$Brush
h2o_ed3092_Activity <- h2o_ed3092$Activity
h2o_ed3092_Travel <- h2o_ed3092$Travel
h2o_ed3092_Curling <- h2o_ed3092$Curling
h2o_ed3092_Asymmetry <- h2o_ed3092$Asymmetry
h2o_ed3092_Stretch <- h2o_ed3092$Stretch
h2o_ed3092_Body <- h2o_ed3092$Body

ed3092_manova <- manova(cbind(h2o_ed3092_Wave,h2o_ed3092_Brush,h2o_ed3092_Activity,h2o_ed3092_Travel,h2o_ed3092_Curling,h2o_ed3092_Asymmetry,h2o_ed3092_Stretch,h2o_ed3092_Body) ~ Age, data = h2o_ed3092)
#summary(ed3092_manova)

#hk104
h2o_hk104_Wave <- h2o_hk104$Wave
h2o_hk104_Brush <- h2o_hk104$Brush
h2o_hk104_Activity <- h2o_hk104$Activity
h2o_hk104_Travel <- h2o_hk104$Travel
h2o_hk104_Curling <- h2o_hk104$Curling
h2o_hk104_Asymmetry <- h2o_hk104$Asymmetry
h2o_hk104_Stretch <- h2o_hk104$Stretch
h2o_hk104_Body <- h2o_hk104$Body

hk104_manova <- manova(cbind(h2o_hk104_Wave,h2o_hk104_Brush,h2o_hk104_Activity,h2o_hk104_Travel,h2o_hk104_Curling,h2o_hk104_Asymmetry,h2o_hk104_Stretch,h2o_hk104_Body) ~ Age, data = h2o_hk104)

#ju775
h2o_ju775_Wave <- h2o_ju775$Wave
h2o_ju775_Brush <- h2o_ju775$Brush
h2o_ju775_Activity <- h2o_ju775$Activity
h2o_ju775_Travel <- h2o_ju775$Travel
h2o_ju775_Curling <- h2o_ju775$Curling
h2o_ju775_Asymmetry <- h2o_ju775$Asymmetry
h2o_ju775_Stretch <- h2o_ju775$Stretch
h2o_ju775_Body <- h2o_ju775$Body

ju775_manova <- manova(cbind(h2o_ju775_Wave,h2o_ju775_Brush,h2o_ju775_Activity,h2o_ju775_Travel,h2o_ju775_Curling,h2o_ju775_Asymmetry,h2o_ju775_Stretch,h2o_ju775_Body) ~ Age, data = h2o_ju775)

#my16
h2o_my16_Wave <- h2o_my16$Wave
h2o_my16_Brush <- h2o_my16$Brush
h2o_my16_Activity <- h2o_my16$Activity
h2o_my16_Travel <- h2o_my16$Travel
h2o_my16_Curling <- h2o_my16$Curling
h2o_my16_Asymmetry <- h2o_my16$Asymmetry
h2o_my16_Stretch <- h2o_my16$Stretch
h2o_my16_Body <- h2o_my16$Body

my16_manova <- manova(cbind(h2o_my16_Wave,h2o_my16_Brush,h2o_my16_Activity,h2o_my16_Travel,h2o_my16_Curling,h2o_my16_Asymmetry,h2o_my16_Stretch,h2o_my16_Body) ~ Age, data = h2o_my16)

#n2
h2o_n2_Wave <- h2o_n2$Wave
h2o_n2_Brush <- h2o_n2$Brush
h2o_n2_Activity <- h2o_n2$Activity
h2o_n2_Travel <- h2o_n2$Travel
h2o_n2_Curling <- h2o_n2$Curling
h2o_n2_Asymmetry <- h2o_n2$Asymmetry
h2o_n2_Stretch <- h2o_n2$Stretch
h2o_n2_Body <- h2o_n2$Body

n2_manova <- manova(cbind(h2o_n2_Wave,h2o_n2_Brush,h2o_n2_Activity,h2o_n2_Travel,h2o_n2_Curling,h2o_n2_Asymmetry,h2o_n2_Stretch,h2o_n2_Body) ~ Age, data = h2o_n2)


sink(paste("GST_CeleST_manova.txt",sep=""), append = FALSE)

writeLines("MANOVA of the Eight Measurements used in LDA\n")

writeLines("########")
writeLines("\nAF16 MANOVA\n")
summary(af16_manova, test = "Pillai")
summary(af16_manova, test = "Wilks")
summary(af16_manova, test = "Hotelling-Lawley")
summary(af16_manova, test = "Roy")

writeLines("########")
writeLines("\nED3092 MANOVA\n")
summary(ed3092_manova, test = "Pillai")
summary(ed3092_manova, test = "Wilks")
summary(ed3092_manova, test = "Hotelling-Lawley")
summary(ed3092_manova, test = "Roy")

writeLines("########")
writeLines("\nHK104 MANOVA\n")
summary(hk104_manova, test = "Pillai")
summary(hk104_manova, test = "Wilks")
summary(hk104_manova, test = "Hotelling-Lawley")
summary(hk104_manova, test = "Roy")

writeLines("########")
writeLines("\nJU775 MANOVA\n")
summary(ju775_manova, test = "Pillai")
summary(ju775_manova, test = "Wilks")
summary(ju775_manova, test = "Hotelling-Lawley")
summary(ju775_manova, test = "Roy")

writeLines("########")
writeLines("\nMY16 MANOVA\n")
summary(my16_manova, test = "Pillai")
summary(my16_manova, test = "Wilks")
summary(my16_manova, test = "Hotelling-Lawley")
summary(my16_manova, test = "Roy")

writeLines("########")
writeLines("\nn2 MANOVA\n")
summary(n2_manova, test = "Pillai")
summary(n2_manova, test = "Wilks")
summary(n2_manova, test = "Hotelling-Lawley")
summary(n2_manova, test = "Roy")

sink()


sink(paste("GST_CeleST_manova_pillai.txt",sep=""), append = FALSE)

writeLines("MANOVA of the Eight Measurements used in LDA\n")

writeLines("########")
writeLines("\nAF16 MANOVA\n")
summary(af16_manova, test = "Pillai")

writeLines("########")
writeLines("\nED3092 MANOVA\n")
summary(ed3092_manova, test = "Pillai")

writeLines("########")
writeLines("\nHK104 MANOVA\n")
summary(hk104_manova, test = "Pillai")

writeLines("########")
writeLines("\nJU775 MANOVA\n")
summary(ju775_manova, test = "Pillai")

writeLines("########")
writeLines("\nMY16 MANOVA\n")
summary(my16_manova, test = "Pillai")

writeLines("########")
writeLines("\nn2 MANOVA\n")
summary(n2_manova, test = "Pillai")

sink()


# plot  -------------------------------------------------------------------

library(ggplot2)
library(ggh4x)
library(introdataviz)
library(dplyr)

celest <- read_excel("GST_CeleST_Output.xlsx")

# Add AgeGroup column
celest <- celest %>%
  mutate(AgeGroup = ifelse(Age %in% c(6, 8), "young", 
                           ifelse(Age %in% c(12, 16), "old", NA)))

celest$AgeGroup <- factor(celest$AgeGroup, levels = c("young", "old"))

celest$StrainAge <- droplevels(with(celest, interaction(Strain,Age, sep="x")))
celest$AgeComp <- droplevels(with(celest, interaction(Age,Compound, sep="x")))

Control <- "CTRL_H2O"
Intervention <- "gold_sodium_thiomalate"

# Generate the plot
celest_plot <- ggplot(
  celest, 
  aes(
    x = factor(AgeGroup), # Age on the x-axis
    y = AdjSwimmingScore, 
    fill = Compound # Split violin based on Compound
  )
) +
  geom_split_violin(scale = "width") + # Create split violin plots
  geom_point(
    aes(color = "black"), 
    pch = 16, 
    size = 0.5, 
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5) # Jitter for individual points
  ) +
  # Add mean lines
  stat_summary(
    fun = "mean", # Compute
    geom = "crossbar", # Draw a horizontal line for the mean
    position = position_dodge(width = 0.9), # Align with the violins
    width = 0.3, # Line width
    color = "black" # Line color
  ) +
  # Add connecting lines for means
  stat_summary(
    fun = "mean",
    geom = "line",
    aes(group = Compound),
    position = position_dodge(width = 0.9),
    color = "black"
  ) +
  scale_fill_manual(
    values =  c("CTRL_H2O" = "#CCCCCC", "gold_sodium_thiomalate" = "gold"),
    name = "Compound" # Legend title for fill
  ) +
  scale_color_manual(
    values =  c("CTRL_H2O" = "#CCCCCC", "gold_sodium_thiomalate" = "gold"),
    name = "Compound" # Legend title for points
  ) +
  scale_y_continuous() +
  ggtitle("GST CeleST: Adjusted Swimming Scores by Age and Compound") +
  xlab("Age") +
  ylab("Adjusted Swimming Score") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels
    legend.position = "right", # Add legend to the right
    legend.title = element_text(size = 14), # Adjust legend title size
    legend.text = element_text(size = 12) # Adjust legend text size
  ) +
  ggh4x::facet_wrap2(
    ~ factor(Strain, levels = c("N2", "JU775", "MY16", "AF16", "ED3092", "HK104")) # Order strains
  )

# Print the plot
print(celest_plot)


# Print the plot
print(celest_plot)


plot(xy)

ggsave("bromo_aging_clock.png", width=8, height=9, xy, device="png")
ggsave("bromo_aging_clock.svg", width=8, height=9, xy, device="svg")




# R version 4.3.3
# car version 3.1-3
# MASS version 7.3-60.0.1
# Bbmisc 1.13
