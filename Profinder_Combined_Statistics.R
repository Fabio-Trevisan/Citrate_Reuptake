library(readr)
library(ggpubr)
library(tidyverse)
library(plyr)
library(dplyr)
library(rstatix)
library(purrr)
library(agricolae)
library(reshape2)

#Profinder (Agilent) software CSIA statistics for AmminoAcids ####
#INCOMPLETE!!!! UN-USEFUL#

#REPLICATE RESULTS USING UNLABELED AS DATA VARIABLE (95% INSTEAD OF 5%) AND CHECK RELIABILITY OF RESULTS
#CHECK CORRELATION BETWEEN [C] & ENRICHMENT% AND RELATIVE ABUNDANCE & ENRICHMENT%
#RUN 2-WAY, 1-WAY ANOVA AND POST HOC TEST ON TREATMENT AND TIME OF SIGNIFICANT MOLECULES

#Read CSV ####
df <- read.csv("DATA enrichment.csv", sep=";",
                  header=T)
df <- df[,-c(2,3)]


#Dataset Uniformation ####
df$Compound <- factor(df$Compound)
df$Time <- factor(df$Time)
df$Labeling <- factor(df$Labeling)
df$Treatment <- factor(df$Treatment)

vector_Compound <- levels(factor(df$Compound))
vector_Treatment <- levels(factor(df$Treatment))
vector_Time <- levels(factor(df$Time))

#USED LABELED VALUES FOR TEST!!! (% OF 13C labelled)
#create Subset according to Compound - Treatment - Time ####
#Compound
Subset <- lapply(vector_Compound, function(i){ 
  i <- subset(df, Compound == i)
})
names(Subset) <- vector_Compound #add names

#Treatment
Subset_2 <- lapply(vector_Compound, function(m){
  lapply(vector_Treatment, function(i){
    subset(Subset[[m]], Treatment == i)
  })
})
names(Subset_2) <- vector_Compound #add names
for(i in vector_Compound) {
  names(Subset_2[[i]]) <- vector_Treatment
} #add names

#Time
Subset_3 <- lapply(vector_Compound, function(m){
  lapply(vector_Treatment, function(i){
    lapply(vector_Time, function(n){
    subset(Subset_2[[m]][[i]], Time == n)
    })
  })
})
names(Subset_3) <- vector_Compound #add names
for(i in vector_Compound) {
  names(Subset_3[[i]]) <- vector_Treatment
} #add names
for(i in vector_Compound) {
  for(n in vector_Treatment) {
    names(Subset_3[[i]][[n]]) <- vector_Time
  }
} #add names



#Summary Table ####
Summary_table_Labeled <- ddply(df, c("Compound", "Time", "Labeling", "Treatment"), summarise,
                               N    = sum(!is.na(Labeled)),
                               mean = mean(Labeled, na.rm=TRUE),
                               sd   = sd(Labeled, na.rm=TRUE),
                               se   = sd / sqrt(N))
write.table(Summary_table_Labeled, file = "Labeled_Summary_table.csv", quote = FALSE, sep = ";")

Enrichment_mean <- dcast(Summary_table_Labeled, Compound + Time + Treatment ~ Labeling, value.var = "mean") #cast (inverse of melt) to reorganice table and be able tu find difference in labeling
Enrichment_se <- dcast(Summary_table_Labeled, Compound + Time + Treatment ~ Labeling, value.var = "se") #cast (inverse of melt) to reorganice table and be able tu find "se" in labeling

Enrichment_df <- Enrichment_mean #combine means and se for both L and U
Enrichment_df$L_se <- Enrichment_se$L
Enrichment_df$U_se <- Enrichment_se$U
Enrichment_df <- Enrichment_df[, c(1, 2, 3, 4, 6, 5, 7)] #reorder columns
Enrichment_df$Enrichment <- Enrichment_df$L - Enrichment_df$U #calculata delta enrichment



#Assumptions ####
## 1. Homogeneity of variances
##Labeling*Time*Treatment
Levene_test_Relative_abundance <- lapply(split(df, df$Compound), function(i){
  levene_test(Relative_abundance ~ Labeling * Time * Treatment, data = i)
})
sink("RelAbund_Levene_test_Homogeneity.csv")
Levene_test_Relative_abundance
sink(NULL)

Levene_test_Unlabeled <- lapply(split(df, df$Compound), function(i){
  levene_test(Unlabeled ~ Labeling * Time * Treatment, data = i)
})
sink("Unlabeled_Levene_test_Homogeneity.csv")
Levene_test_Unlabeled
sink(NULL)

Levene_test_Labeled <- lapply(split(df, df$Compound), function(i){
  levene_test(Labeled ~ Labeling * Time * Treatment, data = i)
})
sink("Labeled_Levene_test_Homogeneity.csv")
Levene_test_Labeled
sink(NULL)

##2. Normality
##Shapiro-Wilk test for all single Treatments
SW_test_RelAbund <- df %>%
  group_by(Labeling, Time, Treatment, Compound) %>%
  shapiro_test(Relative_abundance)
write.table(SW_test_RelAbund, file = "RelAbund_ShapiroWilk_test_Normality.csv", quote = FALSE, sep = ";")

SW_test_Unlabeled <- df %>%
  group_by(Labeling, Time, Treatment, Compound) %>%
  shapiro_test(Unlabeled)
write.table(SW_test_Unlabeled, file = "Unlabeled_ShapiroWilk_test_Normality.csv", quote = FALSE, sep = ";")

SW_test_Labeled <- df %>%
  group_by(Labeling, Time, Treatment, Compound) %>%
  shapiro_test(Labeled)
write.table(SW_test_Labeled, file = "Labeled_ShapiroWilk_test_Normality.csv", quote = FALSE, sep = ";")

##3. Indipendency
#Data are indepent by experimental design!
  
  
  
#3way ANOVA ####
##Multiple Compound Labeled
ThreeWay_Anova <- lapply(split(df, df$Compound), function(i){
  anova(lm(Labeled ~ Labeling * Treatment * Time, data = i))
})
write.table(ThreeWay_Anova, file = "ThreeWay_Anova_.csv", quote = FALSE, sep = ";")

sink("ThreeWay_Anova_2.0_.csv")
ThreeWay_Anova
sink(NULL)

#1way ANOVA ####
#Treatment and Time statistics####
TO BE DONE



#T-Test (wilcox.test) of labeled vs UNlabeled!!! #### 
# Wilcoxon Mann Withney U-test or Wilcoxon Rank sum test for INDIPENDENT data
# not the Wilcoxon sign test for DEPENDENT samples
#greater --> for testing L > U (REAL HYPOTHESIS)
Wilcox_test_greater <- lapply(vector_Compound, function(m){
  lapply(names(Subset_3[[m]]), function(i){ 
    lapply(names(Subset_3[[m]][[i]]), function(n){ 
    wilcox.test(Labeled ~ Labeling, data = Subset_3[[m]][[i]][[n]], alternative = "greater")
    })
  })
})
names(Wilcox_test_greater) <- vector_Compound #add names
for(i in vector_Compound) {
  names(Wilcox_test_greater[[i]]) <- vector_Treatment
} #add names
for(i in vector_Compound) {
  for(n in vector_Treatment) {
    names(Wilcox_test_greater[[i]][[n]]) <- vector_Time
  }
} #add names

#less --> for testing U > L 
Wilcox_test_less <- lapply(vector_Compound, function(m){
  lapply(names(Subset_3[[m]]), function(i){ 
    lapply(names(Subset_3[[m]][[i]]), function(n){ 
      wilcox.test(Labeled ~ Labeling, data = Subset_3[[m]][[i]][[n]], alternative = "less")
    })
  })
})
names(Wilcox_test_less) <- vector_Compound
for(i in vector_Compound) {
  names(Wilcox_test_less[[i]]) <- vector_Treatment
} #add names
for(i in vector_Compound) {
  for(n in vector_Treatment) {
    names(Wilcox_test_less[[i]][[n]]) <- vector_Time
  }
} #add names



#P-Value print in original table ####
##P-Value df preparation ####
Wilcox_test_greater_PValue <- lapply(vector_Compound, function(m){
  lapply(names(Wilcox_test_greater[[m]]), function(i){ 
    lapply(names(Wilcox_test_greater[[m]][[i]]), function(n){
  as.data.frame(Wilcox_test_greater[[m]][[i]][[n]][["p.value"]])
    })
  })
})
names(Wilcox_test_greater_PValue) <- vector_Compound #extract P_Value as list
for(i in vector_Compound) {
  names(Wilcox_test_greater_PValue[[i]]) <- vector_Treatment
} #add names
for(i in vector_Compound) {
  for(n in vector_Treatment) {
    names(Wilcox_test_greater_PValue[[i]][[n]]) <- vector_Time
  }
} #add names

for(i in vector_Compound) {
  for(n in vector_Treatment) {
    for(m in vector_Time) {
      names(Wilcox_test_greater_PValue[[i]][[n]][[m]]) <- "P_Value"
    }
  }
} #rename sub-title (columns) to P-Value



##Transform list in dataframe ####
P_Value <- lapply(vector_Compound, function(m){
  lapply(names(Wilcox_test_greater[[m]]), function(i){ 
    do.call(rbind.data.frame, Wilcox_test_greater_PValue[[m]][[i]])
  })
}) #merge P-Values for Time
names(P_Value) <- vector_Compound #add names
for(i in vector_Compound) {
  names(P_Value[[i]]) <- vector_Treatment
} #add names

P_Value_2 <- lapply(vector_Compound, function(m){
  data.frame(P_Value[[m]])
}) #merge P-Values for Treatment
names(P_Value_2) <- vector_Compound #add names
for(i in vector_Compound) {
  colnames(P_Value_2[[i]]) <- c("C", "Fe", "P")
} #rename sub-title (columns) according to Treatment

for(i in vector_Compound) {
  P_Value_2[[i]]$Time <- rownames(P_Value_2[[i]])
} #move Time from rownames to column

P_Value_3 <- lapply(vector_Compound, function(m){
  melt(P_Value_2[[m]], id = "Time", value.name = "P_Value", variable.name = "Treatment")
}) #melt to bring Treatment from columns to rows
names(P_Value_3) <- vector_Compound #add names

P_Value_df <- do.call(rbind.data.frame, P_Value_3) #merge all sub-dfs
P_Value_df$Compound <- rownames(P_Value_df) #move row nnames to column

P_Value_df$Compound <-gsub("\\.|0|1|2|3|4|5|6|7|8|9","",as.character(P_Value_df$Compound)) #clean column names from ecxessive strings 
P_Value_df <- P_Value_df[, c(4, 1, 2, 3)] #reorder columns
P_Value_df <- P_Value_df[order(P_Value_df$Compound, P_Value_df$Time),] #reorder rows to match Enrichment_df
rownames(P_Value_df) <- NULL #rename Rows



##Combining df ####
Enrichment_df$P_Value <- P_Value_df$P_Value  #Combine df --> add column with P_Values to original table

##P Value filtering ####
z <- 0.05 #filtering ratio
Filtered <- Enrichment_df[Enrichment_df$`P_Value`<z,] #P_Value filtering >z

#save
write.csv(Filtered, file = "AA_Significant_Enrichment.csv", row.names=FALSE)





#summary other variables () ####
Summary_table_Relative_abundance <- ddply(df, c("Compound", "Time", "Labeling", "Treatment"), summarise,
                                          N    = sum(!is.na(Relative_abundance)),
                                          mean = mean(Relative_abundance, na.rm=TRUE),
                                          sd   = sd(Relative_abundance, na.rm=TRUE),
                                          se   = sd / sqrt(N))
write.table(Summary_table_Relative_abundance, file = "RelAbund_Summary_table.csv", quote = FALSE, sep = ";")

Summary_table_Unlabeled <- ddply(df, c("Compound", "Time", "Labeling", "Treatment"), summarise,
                                 N    = sum(!is.na(Unlabeled)),
                                 mean = mean(Unlabeled, na.rm=TRUE),
                                 sd   = sd(Unlabeled, na.rm=TRUE),
                                 se   = sd / sqrt(N))
write.table(Summary_table_Unlabeled, file = "Unlabeled_Summary_table.csv", quote = FALSE, sep = ";")
