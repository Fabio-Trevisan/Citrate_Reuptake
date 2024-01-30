library(readr)
library(ggpubr)
library(tidyverse)
library(plyr)
library(dplyr)
library(rstatix)
library(purrr)
library(agricolae)
library(reshape2)


#CHECK CORRELATION BETWEEN RELATIVE ABUNDANCE & ENRICHMENT%

#OA
#Read CSV ####
table <- read.csv("20230711 3-NPH acids from in-house script_IsoCor_res.tsv", sep="\t", header=T)
Class <- "OA"

#AA
#Read CSV ####
table <- read.csv("20230714_AccQ-Tag_AA_peak areas_inhouse_script 2.0_IsoCor_res.tsv", sep="\t", header=T)
Class <- "AA"

#Dataset preparation ####
#cleaning
df <- table[,-c(3,5,7:9)] #remove un-useful columns
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

#Dataset Uniformation ####
df$metabolite <- factor(df$metabolite)
df$Time <- factor(df$Time)
df$Labeling <- factor(df$Labeling)
df$Treatment <- factor(df$Treatment)

vector_metabolite <- levels(factor(df$metabolite))
vector_Treatment <- levels(factor(df$Treatment))
vector_Time <- levels(factor(df$Time))

#Replace 0
df["mean_enrichment"][df["mean_enrichment"]==0] <- NA #raplace meanenrichment 0 values with 1/10 ot 1/100 of min values
x <- min(df$mean_enrichment, na.rm = T)/10 #calculate min values/10
df[is.na(df)] <- x #raplace meanenrichment 0 values with 1/10 of min values

#create Subset according to metabolite - Treatment - Time ####
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



#Correlation statistic test ####
res <- lapply(vector_metabolite, function(m){
  lapply(names(Subset_3[[m]]), function(i){ 
    lapply(names(Subset_3[[m]][[i]]), function(n){ 
      cor.test(Subset_3[[m]][[i]][[n]][["area"]], Subset_3[[m]][[i]][[n]][["mean_enrichment"]], method = "spearman")
    })
  })
})
names(res) <- vector_metabolite #add names
for(i in vector_metabolite) {
  names(res[[i]]) <- vector_Treatment
} #add names
for(i in vector_metabolite) {
  for(n in vector_Treatment) {
    names(res[[i]][[n]]) <- vector_Time
  }
} #add names

##print all ####
sink(paste("InHouse_", Class, "_Corr.Table.csv", sep=""))
res
sink(NULL)


##print only P-Value and Rho####
#Rho
Rho <- lapply(vector_metabolite, function(m){
  lapply(names(res[[m]]), function(i){ 
    lapply(names(res[[m]][[i]]), function(n){ 
      as.data.frame(res[[m]][[i]][[n]][["estimate"]])
    })
  })
}) #merge rho estimate
names(Rho) <- vector_metabolite #add names
for(i in vector_metabolite) {
  names(Rho[[i]]) <- vector_Treatment
} #add names
for(i in vector_metabolite) {
  for(n in vector_Treatment) {
    names(Rho[[i]][[n]]) <- vector_Time
  }
} #add names

Rho_1 <- lapply(vector_metabolite, function(m){
  lapply(names(Rho[[m]]), function(i){ 
    do.call(rbind.data.frame, Rho[[m]][[i]])
  })
}) #merge P Values for time
names(Rho_1) <- vector_metabolite #add names
for(i in vector_metabolite) {
  names(Rho_1[[i]]) <- vector_Treatment
} #add names

Rho_2 <- lapply(vector_metabolite, function(m){
  data.frame(Rho_1[[m]])
}) #merge P-Values for Treatment
names(Rho_2) <- vector_metabolite #add names
for(i in vector_metabolite) {
  colnames(Rho_2[[i]]) <- c("C", "Fe", "P")
} #rename sub-title (columns) according to Treatment
for(i in vector_metabolite) {
  Rho_2[[i]]$Time <- rownames(Rho_2[[i]])
} #move Time from rownames to column

Rho_3 <- lapply(vector_metabolite, function(m){
  melt(Rho_2[[m]], id = "Time", value.name = "Rho", variable.name = "Treatment")
}) #melt to bring Treatment from columns to rows
names(Rho_3) <- vector_metabolite #add names

Rho_df <- do.call(rbind.data.frame, Rho_3) #merge all sub-dfs
Rho_df$metabolite <- rownames(Rho_df) #move row nnames to column

Rho_df$metabolite <-gsub("\\.|0|1|2|3|4|5|6|7|8|9","",as.character(Rho_df$metabolite)) #clean column names from ecxessive strings 
Rho_df <- Rho_df[, c(4, 1, 2, 3)] #reorder columns
Rho_df$Time <- as.integer(Rho_df$Time)
Rho_df <- Rho_df[order(Rho_df$metabolite, Rho_df$Time),] #re-order rows to match Enrichment_df
rownames(Rho_df) <- NULL #rename Rows



#P-Value
P_Value <- lapply(vector_metabolite, function(m){
  lapply(names(res[[m]]), function(i){ 
    lapply(names(res[[m]][[i]]), function(n){ 
      as.data.frame(res[[m]][[i]][[n]][["p.value"]])
    })
  })
}) #extract P Values and transform in df 
names(P_Value) <- vector_metabolite #add names
for(i in vector_metabolite) {
  names(P_Value[[i]]) <- vector_Treatment
} #add names
for(i in vector_metabolite) {
  for(n in vector_Treatment) {
    names(P_Value[[i]][[n]]) <- vector_Time
  }
} #add names

P_Value_1 <- lapply(vector_metabolite, function(m){
  lapply(names(P_Value[[m]]), function(i){ 
    do.call(rbind.data.frame, P_Value[[m]][[i]])
  })
}) #merge P Values for time
names(P_Value_1) <- vector_metabolite #add names
for(i in vector_metabolite) {
  names(P_Value_1[[i]]) <- vector_Treatment
} #add names

P_Value_2 <- lapply(vector_metabolite, function(m){
  data.frame(P_Value_1[[m]])
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



##Combining df ####
Spearman <- P_Value_df 
Spearman$Rho <- Rho_df$Rho #Combine df --> add column with Rho to P-Value table
Spearman <- na.omit(Spearman)

##P Value filtering ####
z <- 0.05 #filtering ratio
Filtered <- Spearman[Spearman$`P_Value`>z,] #P_Value filtering >z

#save
write.table(Spearman, file = paste("InHouse_", Class, "_Corr.Table_2.0.csv", sep=""), row.names=FALSE, sep = ";")
write.table(Filtered, file = paste("InHouse_", Class, "_Corr.Table_2.0_(NOT).csv", sep=""), row.names=FALSE, sep = ";")



# Correlation plots####
#Time as colour + Treatment as FacetWrap 
p1 <- ggscatter(df, x = "mean_enrichment", y = "area", 
          color = "Time",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Enrichment %", ylab = "Relative_abundance") +
  scale_color_manual(values=c("grey77","skyblue3", "blue", "darkblue")) +
  scale_fill_manual(values=c("grey77","skyblue3", "blue", "darkblue")) +
  facet_wrap(~metabolite + Treatment, scales="free")

ggsave(filename = paste("InHouse_", Class, "_Corr-matrix_(Time).pdf", sep=""), plot = p1, dpi = 600, units = "cm", width = 150, height = 110, scale = 0.5)

#Treatment as colour + Time as FacetWrap 
p2 <- ggscatter(df, x = "mean_enrichment", y = "area", 
               color = "Treatment",
               add = "reg.line", conf.int = TRUE, 
               cor.coef = TRUE, cor.method = "spearman",
               xlab = "Enrichment %", ylab = "Relative_abundance") +
  scale_color_manual(values=c("grey77","orange", "blue")) +
  scale_fill_manual(values=c("grey77","orange", "blue")) +
  facet_wrap(~metabolite + Time, scales="free")

ggsave(filename = paste("InHouse_", Class, "_Corr-matrix_(Treatment).pdf", sep=""), plot = p2, dpi = 600, units = "cm", width = 150, height = 110, scale = 0.5)

#Metabolites as FacetWrap 
p3 <- ggscatter(df, x = "mean_enrichment", y = "area", 
               add = "reg.line", conf.int = TRUE, 
               cor.coef = TRUE, cor.method = "spearman",
               xlab = "Enrichment %", ylab = "Relative_abundance") +
  facet_wrap(~metabolite, scales="free")

ggsave(filename = paste("InHouse_", Class, "_Corr-matrix.pdf", sep=""), plot = p3, dpi = 600, units = "cm", width = 150, height = 110, scale = 0.5)

#All together
p4 <- ggscatter(df, x = "mean_enrichment", y = "area", 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "spearman",
                xlab = "Enrichment %", ylab = "Relative_abundance")

ggsave(filename = paste("InHouse_", Class, "_Correlation.pdf", sep=""), plot = p4, dpi = 600, units = "cm", width = 75, height = 55, scale = 0.3)
