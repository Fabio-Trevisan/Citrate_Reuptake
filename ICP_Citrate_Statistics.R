library(readr)
library(ggpubr)
library(tidyverse)
library(plyr)
library(dplyr)
library(rstatix)
library(purrr)
library(agricolae)
library(reshape2)


#ICP statistics ####
#Read CSV ####
df <- read.csv("DATA_Citrate_ICP.csv", sep=";",
                  header=T)

df2 <- melt(data = table, id.vars = c("Tissue", "Treatment", "Concentration"), 
               variable.name = "Element", 
               value.name = "ppm")


Summary_table <- ddply(df2, c("Tissue", "Treatment", "Element"), summarise,
                       N    = sum(!is.na(ppm)),
                       mean = mean(ppm, na.rm=TRUE),
                       sd   = sd(ppm, na.rm=TRUE),
                       se   = sd / sqrt(N))

ICP_table <- Summary_table[,-6]
ICP_table$mean <- round(ICP_table$mean, 3)
ICP_table$se <- round(ICP_table$se, 3)
ICP_table$Concentration <- paste(ICP_table$mean, ICP_table$se, sep = " +/- ")
ICP_table <- ICP_table[,-c(5,6)]

ICP_table_2 <- dcast(ICP_table, Tissue+Treatment+N~Element, value.var = "Concentration") 


write.table(ICP_table_2, file = "ICP_table_Citrate.csv", quote = FALSE, sep = ";")



#Assumptions ####
## 1. Homogeneity of variances
##Treatment*Tissue
Levene_test2 <- df2 %>%
  group_by(Element) %>%
  levene_test(ppm ~ Treatment * Tissue)

##2. Normality
##Shapiro-Wilk test for all single treatments
SW_test <- df2 %>%
  group_by(Tissue, Element, Treatment) %>%
  shapiro_test(ppm)
View(SW_test)
write.table(SW_test, file = "ICP_ShapiroWilk_Citrate.csv", quote = FALSE, sep = ";")

##3. Indipendency
#Data are indepent by experimental design!


#1way ANOVA ####
###create Subsets according to Tissue 
vector_Tissue <- c("R", "S")

Subsets <- lapply(vector_Tissue, function(i){ 
  i <- subset(df2, Tissue == i)
})

names(Subsets) <- vector_Tissue

##Treatment.for tukey
OneWay_Anova_Tr <- lapply(vector_Tissue, function(m){
  lapply(split(Subsets[[m]], Subsets[[m]][["Element"]]), function(i){ 
    aov(ppm ~ Treatment, data = i)
  })
})
names(OneWay_Anova_Tr) <- vector_Tissue


##Treatment.for print
OneWay_Anova_Tr2 <- lapply(vector_Tissue, function(m){
  lapply(split(Subsets[[m]], Subsets[[m]][["Element"]]), function(i){ 
    anova(lm(ppm ~ Treatment, data = i))
  })
})
names(OneWay_Anova_Tr2) <- vector_Tissue


##OneWayAnova save
sink("ICP_OneWayAnova_Citrate_Tr.csv")
OneWay_Anova_Tr2 
sink(NULL)


#Tukey as post hoc test ####
##Treatment
HSD_Tr <- lapply(vector_Tissue, function(m){
  lapply(names(OneWay_Anova_Tr[[m]]), function(i){ 
    HSD.test(OneWay_Anova_Tr[[m]][[i]], "Treatment")
  })
})
names(HSD_Tr) <- vector_Tissue
for(i in vector_Tissue) {
  list <- names(OneWay_Anova_Tr[[i]]) 
  names(HSD_Tr[[i]]) <- list
}


##HSD_test save
##Treatment
HSD_Tr_groups <- lapply(vector_Tissue, function(i){
  lapply(names(OneWay_Anova_Tr[[i]]), function(m){
    as.data.frame(HSD_Tr[[i]][[m]][["groups"]])
  })
})
names(HSD_Tr_groups) <- vector_Tissue
for(i in vector_Tissue) {
  list <- names(OneWay_Anova_Tr[[i]]) 
  names(HSD_Tr_groups[[i]]) <- list
}
sink("ICP_HSD_Citrate_Tr.csv")
HSD_Tr_groups 
sink(NULL)

