library(readr)
library(ggpubr)
library(ggplot2)
library(tidyverse)
library(plyr)
library(dplyr)
library(rstatix)
library(purrr)
library(agricolae)
library(reshape2)
library(scales)


#OA
# Read CSV ####
table <- read.csv("20230711 3-NPH acids from in-house script_IsoCor_res.tsv", sep="\t", header=T)
Class <- "OA"

#AA
# Read CSV ####
table <- read.csv("20230714_AccQ-Tag_AA_peak areas_inhouse_script 2.0_IsoCor_res.tsv", sep="\t", header=T)
Class <- "AA"

# Dataset preparation ####
#cleaning
df <- table[,-c(3,5:9)] #remove un-useful columns
df <- df[df$isotopologue<=0,] #filter only isotopologue M+0
df <- df[!str_detect(df$sample, "blank"),] #remove blank samples
df <- df[!str_detect(df$sample, "std"),] #remove standards
df <- df[!str_detect(df$sample, "QC"),] #remove QCs
row.names(df) <- NULL
df <- df[,-3] #remove un-useful columns

#Time-Treatment-Labeling 
Treatment_factors <- read.csv(paste("Treatment_factors_", Class, ".csv", sep=""), sep=";", header=T) #load file with treatment_factors
df$Treatment <- Treatment_factors$Treatment
df$Time <- Treatment_factors$Time
df$Labeling <- Treatment_factors$Labeling #add them to df

#Uniformation and vectors creation
df$metabolite <- factor(df$metabolite)
df$Treatment <- factor(df$Treatment)
df$Time <- factor(df$Time)
df$Labeling <- factor(df$Labeling) #convert character and numerbers to factors

vector_metabolite <- levels(factor(df$metabolite))
vector_Treatment <- levels(factor(df$Treatment))
vector_Time <- levels(factor(df$Time))
vector_Labeling <- levels(factor(df$Labeling)) #write vectors for subseting

#Replace 0 with random values of 1/10 +/- 10% of min value
df["mean_enrichment"][df["mean_enrichment"]==0] <- NA #raplace meanenrichment 0 values with 1/10 of min values
x <- min(df$mean_enrichment, na.rm = T)/10 #calculate min values/10
df[is.na(df)] <-  runif(sum(is.na(df)), min = x-(x/10), max = x+(x/10))



##    Subset according to metabolite - Treatment - Time ####
#metabolite
Subset <- lapply(vector_metabolite, function(i){ 
  i <- subset(df, metabolite == i)
})
names(Subset) <- vector_metabolite #add names

#Treatment
Subset_2 <- lapply(vector_metabolite, function(m){
  lapply(vector_Treatment, function(i){
    subset(Subset[[m]], Treatment == i)
  })
})
names(Subset_2) <- vector_metabolite #add names
for(i in vector_metabolite) {
  names(Subset_2[[i]]) <- vector_Treatment
} #add names

#Time
Subset_3 <- lapply(vector_metabolite, function(m){
  lapply(vector_Treatment, function(i){
    lapply(vector_Time, function(n){
    subset(Subset_2[[m]][[i]], Time == n)
    })
  })
})
names(Subset_3) <- vector_metabolite #add names
for(i in vector_metabolite) {
  names(Subset_3[[i]]) <- vector_Treatment
} #add names
for(i in vector_metabolite) {
  for(n in vector_Treatment) {
    names(Subset_3[[i]][[n]]) <- vector_Time
  }
} #add names



# Assumptions ####
## 1. Homogeneity of variances
##Labeling*Time*Treatment
Levene_test <- lapply(split(df, df$metabolite), function(i){
  levene_test(mean_enrichment ~ Labeling * Time * Treatment, data = i)
})
sink(paste("InHouse_", Class, "_Levene_test_Homogeneity.csv", sep=""))
Levene_test
sink(NULL)

##2. Normality
##Shapiro-Wilk test for all single Treatments 
SW_test_Labeled <- df %>%
  group_by(Labeling, Time, Treatment, metabolite) %>%
  shapiro_test(mean_enrichment)
write.table(SW_test_Labeled, file = paste("InHouse_", Class, "_ShapiroWilk_test_Normality.csv", sep=""), quote = FALSE, sep = ";")

##3. Indipendency
#Data are indepent by experimental design!



# 3way ANOVA ####
##for each Metabolite test treatment, time and LvsU
ThreeWay_Anova <- lapply(split(df, df$metabolite), function(i){
  anova(lm(mean_enrichment ~ Labeling * Treatment * Time, data = i))
})
write.table(ThreeWay_Anova, file = paste("InHouse_", Class, "_ThreeWay_Anova.csv", sep=""), quote = FALSE, sep = ";")

sink(paste("InHouse_", Class, "_ThreeWay_Anova_2.0.csv", sep=""))
ThreeWay_Anova
sink(NULL)



# 2way ANOVA on Labeled samples ####
Subset_L <- lapply(vector_metabolite, function(m){
  Subset[[m]] %>%
    filter(str_detect(Labeling, "L"))
}) # remove Unlabeled values since stat will be done comparing labeled values/samples
names(Subset_L) <- vector_metabolite #add names

TwoWay_Anova <- lapply(vector_metabolite, function(m){
    anova(lm(mean_enrichment ~ Treatment * Time, data = Subset_L[[m]]))
})
names(TwoWay_Anova) <- vector_metabolite #add names

write.table(TwoWay_Anova, file = paste("InHouse_", Class, "_TwoWay_Anova.csv", sep=""), quote = FALSE, sep = ";")

sink(paste("InHouse_", Class, "_TwoWay_Anova_2.0.csv", sep=""))
TwoWay_Anova
sink(NULL)



# Kruskal Wallis on Labeled samples ####
# for Treatment and Time comparisons
## The post hoc nonparametrics tests (kruskal) are using the criterium Fisher's least significant difference (LSD)
##Treatment
KW_Tr <- lapply(vector_metabolite, function(m){
  lapply(split(Subset_L[[m]], Subset_L[[m]][["Time"]]), function(i){ 
    kruskal(i$mean_enrichment, i$Treatment, p.adj = "fdr")
  })
})
names(KW_Tr) <- vector_metabolite

##Time
KW_Ti <- lapply(vector_metabolite, function(m){
  lapply(split(Subset_L[[m]], Subset_L[[m]][["Treatment"]]), function(i){ 
    kruskal(i$mean_enrichment, i$Time, p.adj = "fdr")
  })
})
names(KW_Ti) <- vector_metabolite

## Save
##Treatment
KW_Tr_groups <- lapply(vector_metabolite, function(i){
  lapply(names(KW_Tr[[i]]), function(m){
    as.data.frame(KW_Tr[[i]][[m]][["groups"]])
  })
})
names(KW_Tr_groups) <- vector_metabolite
for(i in vector_metabolite) {
  list <- names(KW_Tr[[i]]) 
  names(KW_Tr_groups[[i]]) <- list
}
sink(paste("InHouse_", Class, "_KW_Tr.csv", sep=""))
KW_Tr_groups 
sink(NULL)

##Time
KW_Ti_groups <- lapply(vector_metabolite, function(i){
  lapply(names(KW_Ti[[i]]), function(m){
    as.data.frame(KW_Ti[[i]][[m]][["groups"]])
  })
})
names(KW_Ti_groups) <- vector_metabolite
for(i in vector_metabolite) {
  list <- names(KW_Ti[[i]]) 
  names(KW_Ti_groups[[i]]) <- list
}
sink(paste("InHouse_", Class, "_KW_Ti.csv", sep=""))
KW_Ti_groups 
sink(NULL)



# T-Test (wilcox.test) UvsL#### 
# Wilcoxon Mann Withney U-test or Wilcoxon Rank sum test for INDIPENDENT data
# not the Wilcoxon sign test for DEPENDENT samples
#greater --> for testing L > U (REAL HYPOTHESIS)
Wilcox_test_greater <- lapply(vector_metabolite, function(m){
  lapply(names(Subset_3[[m]]), function(i){ 
    lapply(names(Subset_3[[m]][[i]]), function(n){ 
    wilcox.test(mean_enrichment ~ Labeling, data = Subset_3[[m]][[i]][[n]], alternative = "greater")
    })
  })
})
names(Wilcox_test_greater) <- vector_metabolite #add names
for(i in vector_metabolite) {
  names(Wilcox_test_greater[[i]]) <- vector_Treatment
} #add names
for(i in vector_metabolite) {
  for(n in vector_Treatment) {
    names(Wilcox_test_greater[[i]][[n]]) <- vector_Time
  }
} #add names

#less --> for testing U > L 
Wilcox_test_smaller <- lapply(vector_metabolite, function(m){
  lapply(names(Subset_3[[m]]), function(i){ 
    lapply(names(Subset_3[[m]][[i]]), function(n){ 
      wilcox.test(mean_enrichment ~ Labeling, data = Subset_3[[m]][[i]][[n]], alternative = "less")
    })
  })
})
names(Wilcox_test_smaller) <- vector_metabolite
for(i in vector_metabolite) {
  names(Wilcox_test_smaller[[i]]) <- vector_Treatment
} #add names
for(i in vector_metabolite) {
  for(n in vector_Treatment) {
    names(Wilcox_test_smaller[[i]][[n]]) <- vector_Time
  }
} #add names



# P-Value print in original table ####
##  P-Value df preparation ####
Wilcox_test_greater_PValue <- lapply(vector_metabolite, function(m){
  lapply(names(Wilcox_test_greater[[m]]), function(i){ 
    lapply(names(Wilcox_test_greater[[m]][[i]]), function(n){
  as.data.frame(Wilcox_test_greater[[m]][[i]][[n]][["p.value"]])
    })
  })
})
names(Wilcox_test_greater_PValue) <- vector_metabolite #extract P_Value as list
for(i in vector_metabolite) {
  names(Wilcox_test_greater_PValue[[i]]) <- vector_Treatment
} #add names
for(i in vector_metabolite) {
  for(n in vector_Treatment) {
    names(Wilcox_test_greater_PValue[[i]][[n]]) <- vector_Time
  }
} #add names

for(i in vector_metabolite) {
  for(n in vector_Treatment) {
    for(m in vector_Time) {
      names(Wilcox_test_greater_PValue[[i]][[n]][[m]]) <- "P_Value"
    }
  }
} #rename sub-title (columns) to P-Value



##  Transform list in dataframe ####
P_Value <- lapply(vector_metabolite, function(m){
  lapply(names(Wilcox_test_greater[[m]]), function(i){ 
    do.call(rbind.data.frame, Wilcox_test_greater_PValue[[m]][[i]])
  })
}) #merge P-Values for Time
names(P_Value) <- vector_metabolite #add names
for(i in vector_metabolite) {
  names(P_Value[[i]]) <- vector_Treatment
} #add names

P_Value_2 <- lapply(vector_metabolite, function(m){
  data.frame(P_Value[[m]])
}) #merge P-Values for Treatment
names(P_Value_2) <- vector_metabolite #add names
for(i in vector_metabolite) {
  colnames(P_Value_2[[i]]) <- c("C", "Fe", "P")
} #rename sub-title (columns) according to Treatment

for(i in vector_metabolite) {
  P_Value_2[[i]]$Time <- rownames(P_Value_2[[i]])
} #move Time from rownames to column

P_Value_3 <- lapply(vector_metabolite, function(m){
  melt(P_Value_2[[m]], id = "Time", value.name = "P_Value", variable.name = "Treatment")
}) #melt to bring Treatment from columns to rows
names(P_Value_3) <- vector_metabolite #add names

P_Value_df <- do.call(rbind.data.frame, P_Value_3) #merge all sub-dfs
P_Value_df$metabolite <- rownames(P_Value_df) #move row nnames to column

P_Value_df$metabolite <-gsub("\\.|0|1|2|3|4|5|6|7|8|9","",as.character(P_Value_df$metabolite)) #clean column names from ecxessive strings 
P_Value_df <- P_Value_df[, c(4, 1, 2, 3)] #reorder columns
P_Value_df$Time <- as.integer(P_Value_df$Time)
P_Value_df <- P_Value_df[order(P_Value_df$metabolite, P_Value_df$Time),] #re-order rows to match Enrichment_df
rownames(P_Value_df) <- NULL #rename Rows



# Summary Table & Mean Enrichment ####
Summary_table <- ddply(df, c("metabolite", "Time", "Labeling", "Treatment"), summarise,
                       N    = sum(!is.na(mean_enrichment)),
                       mean = mean(mean_enrichment, na.rm=TRUE),
                       sd   = sd(mean_enrichment, na.rm=TRUE),
                       se   = sd / sqrt(N))
write.table(Summary_table, file = paste("InHouse_", Class, "_Summary_table.csv", sep=""), 
            quote = FALSE, sep = ";")


Enrichment_mean <- dcast(Summary_table, metabolite + Time + Treatment ~ Labeling, value.var = "mean") #cast (inverse of melt) to reorganize table and be able to find difference in labeling
Enrichment_se <- dcast(Summary_table, metabolite + Time + Treatment ~ Labeling, value.var = "se") #cast (inverse of melt) to reorganize table and be able to find "se" in labeling

Enrichment_df <- Enrichment_mean #combine means and se for both L and U
Enrichment_df$L_se <- Enrichment_se$L
Enrichment_df$U_se <- Enrichment_se$U
Enrichment_df <- Enrichment_df[, c(1, 2, 3, 4, 6, 5, 7)] #reorder columns
Enrichment_df$Enrichment <- Enrichment_df$L - Enrichment_df$U #calculate delta enrichment



##  Combining df ####
Enrichment_df$P_Value <- P_Value_df$P_Value  #Combine df --> add column with P_Values to original table

##  P Value filtering ####
z <- 0.05 #filtering ratio
Filtered <- Enrichment_df[Enrichment_df$`P_Value`<z,] #P_Value filtering >z

##  Save ####
write.table(Enrichment_df, file = paste("InHouse_", Class, "_Enrichment.csv", sep=""), row.names=FALSE, sep = ";")
write.table(Filtered, file = paste("InHouse_", Class, "_Significant_Enrichment.csv", sep=""), row.names=FALSE, sep = ";")


