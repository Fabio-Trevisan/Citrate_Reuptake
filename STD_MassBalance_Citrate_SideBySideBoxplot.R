library(readr)
library(ggpubr)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(plyr) 
library(dplyr)
library(readxl)
library(scales)


table <- read.csv("DATA_STD_MassBalance_Citrate.csv", sep=";",
                  header=T)

table <- table %>% drop_na(Value)


#SideBySide Boxlot ####
my_order<- c("C", "Fe", "P")
f1 <- ggplot(table, aes(x = factor(Time), y = Value, fill = Treatment)) +  
  stat_boxplot(geom="errorbar") +
  geom_boxplot() +
  theme_bw() + 
  scale_fill_manual(values=c("grey77","darkorange2", "skyblue3")) +
  scale_color_manual(values=c("grey77","darkorange2", "skyblue3")) +
  scale_y_continuous(labels = scientific)


f2 <- f1 + facet_wrap(~Species_Tissue, scales="free", ncol = 2)+
  ylab("Mass Balance (mg/g Citrate/DW)") + 
  xlab("Time (Days)")
f2
ggsave(filename = "STD_MassBalance_Citrate_SideBySide_Boxplot.pdf", plot = last_plot(), dpi = 600, units = "cm", width = 70, height = 60, scale = 0.5)


#SideBySide Boxlot + trendline ####
trendline <- geom_smooth(aes(group = Treatment, color = Treatment, fill = Treatment), method=lm, alpha = 0.1, linetype="dashed") 

f1 <- ggplot(table, aes(x = factor(Time), y = Value, fill = Treatment)) +  
  trendline +
  stat_boxplot(geom="errorbar") +
  geom_boxplot() +
  theme_bw() + 
  scale_fill_manual(values=c("grey77","darkorange2", "skyblue3")) +
  scale_color_manual(values=c("grey77","darkorange2", "skyblue3")) +
  scale_y_continuous(labels = scientific) +
  stat_regline_equation()

f2 <- f1 + facet_wrap(~Species_Tissue, scales="free", ncol = 2)+
  ylab("Mass Balance (mg/g Citrate/DW)") + 
  xlab("Time (Days)")
f2
ggsave(filename = "STD_MassBalance_Citrate_SideBySide_Boxplot_2.pdf", plot = last_plot(), dpi = 600, units = "cm", width = 70, height = 60, scale = 0.5)


#Regression line statistics ####
##Assumptions  ####
###1. Indipendence of observation (no autocorrelation between time and treatment)
cor(table$Time, table$Treatment)

###2. Normality
hist(table$Value)

###3. Linearity
plot(Value ~ Time, data=table)
plot(Value ~ Treatment, data=table)

###4. Homoscedasticity


##Model fitting and statistics ####
vector_Species_Tissue <- c("50microM_TR",
                           "50microM_TS",
                           "500microM_TR",
                           "500microM_TS")
vector_Treatment <- c("C",
                      "P",
                      "Fe")

###create Subsets according to Species_Tissue ####
Subsets <- lapply(vector_Species_Tissue, function(i){ 
  i <- subset(table, Species_Tissue == i)
})
names(Subsets) <- vector_Species_Tissue


###Regression model development  ####
Regression_line <- lapply(vector_Species_Tissue, function(m){
  lapply(split(Subsets[[m]], Subsets[[m]][["Treatment"]]), function(i){ 
    lm(Value ~ Time, data = i)
  })
})
names(Regression_line) <- vector_Species_Tissue


###Extract statistics/data from model ####
TR50microM <- lapply(vector_Treatment, function(m){
  summary(Regression_line[["50microM_TR"]][[m]])
})
names(TR50microM) <- vector_Treatment
TS50microM <- lapply(vector_Treatment, function(m){
  summary(Regression_line[["50microM_TS"]][[m]])
})
names(TS50microM) <- vector_Treatment
TR500microM <- lapply(vector_Treatment, function(m){
  summary(Regression_line[["500microM_TR"]][[m]])
})
names(TR500microM) <- vector_Treatment
TS500microM <- lapply(vector_Treatment, function(m){
  summary(Regression_line[["500microM_TS"]][[m]])
})
names(TS500microM) <- vector_Treatment

###Save/print results of linear model significance ####
sink("STD_MassBalance_Citrate_Regression_line.csv")
"TR50microM"
TR50microM
"TS50microM"
TS50microM
"TR500microM"
TR500microM
"TS500microM"
TS500microM
sink(NULL)


###Regression_model_comparison  ####
Regression_model_comparison <- lapply(split(table, table$Species_Tissue), function(i){ 
  anova(lm(Value ~ Time*Treatment, data = i))
})

sink("STD_MassBalance_Citrate_Regression_model_Comparison.csv")
Regression_model_comparison
sink(NULL)

