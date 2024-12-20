library(readr)
library(ggpubr)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(plyr) 
library(dplyr)
library(readxl)
library(rstatix)
library(purrr)
library(agricolae)
library(caret)
library(tidyr)


table <- read.csv("DATA_WinRhizor.csv", sep=";",
                  header=T)

table2 <- melt(table, value.name = "Value", id = c("Treatment"), variable.name = "Parameter" )

#Sumamry table####
Summary_table_Winrhizo <- ddply(table2, c("Parameter", "Treatment"), summarise,
                            N    = sum(!is.na(Value)),
                            mean = mean(Value, na.rm=TRUE),
                            sd   = sd(Value, na.rm=TRUE),
                            se   = sd / sqrt(N))
write.table(Summary_table_Winrhizo, file = "WinRhizor_summary_table.csv", quote = FALSE, sep = ";")


#Boxlot ####
my_order<- c("C", "Fe", "P")
b1 <- ggplot(table2, aes(x = factor(Treatment, level = my_order), y = Value, fill = factor(Treatment, level = my_order))) +  
  stat_boxplot(geom="errorbar") +
  geom_boxplot() +
  theme_bw() + 
  scale_fill_manual(values=c("grey77","darkorange2", "skyblue3"), name = "Treatment")
b1

b2 <- b1 + facet_wrap(~Parameter, scales="free", ncol = 3)+
  xlab("Treatment")
b2

ggsave(filename = "WinRhizor_Parameters_Boxplot.pdf", plot = last_plot(), dpi = 600, units = "cm", width = 70, height = 60, scale = 0.5)


##statistics ####
colnames(table2) <- c("Treatment", "Parameter", "Value") 
name <- "WinRhizor_Parameters"
#Assumptions 
## 1. Homogeneity of variances
##Treatment*Time
table2$Treatment <- factor(table2$Treatment)
table2$Parameter <- factor(table2$Parameter)

L_test <- table2 %>%
  group_by(Parameter) %>%
  levene_test(Value ~ Treatment)
View(L_test)
write.table(L_test, file = "WinRhizor_Parameters_Levene_test_results.csv", quote = FALSE, sep = ";")


##2. Normality
##Shapiro-Wilk test for all single treatments
SW_test <- table2 %>%
  group_by(Treatment, Parameter) %>%
  shapiro_test(Value)
View(SW_test)
write.table(SW_test, file = "WinRhizor_Parameters_ShapiroWilk_test_results.csv", quote = FALSE, sep = ";")

##3. Indipendency
Data are indepent by experimental design!
  
##anova for tukey
OneWay_Anova_Boxplot <- lapply(split(table2, table2[["Parameter"]]), function(i){ 
  aov(Value ~ Treatment, data = i)
})

##anova for print
OneWay_Anova_Boxplot2 <- lapply(split(table2, table2[["Parameter"]]), function(i){ 
  anova(lm(Value ~ Treatment, data = i))
})

sink(paste(name, "OneWay_Anova_Boxplot.csv", sep = "_"))
OneWay_Anova_Boxplot2
sink(NULL)

#Tukey
##HSD complete
HSD_Boxplot <- lapply(names(OneWay_Anova_Boxplot), function(i){ 
  HSD.test(OneWay_Anova_Boxplot[[i]], "Treatment")
})
names(HSD_Boxplot) <- names(OneWay_Anova_Boxplot)

##HSD groups only
HSD_Boxplot_groups <- lapply(names(OneWay_Anova_Boxplot), function(i){
  as.data.frame(HSD_Boxplot[[i]][["groups"]])
})
names(HSD_Boxplot_groups) <- names(OneWay_Anova_Boxplot)

sink(paste(name, "HSD_Boxplot.csv", sep = "_"))
HSD_Boxplot_groups
sink(NULL)


#SideBySide Boxlot ####
process1 <- preProcess(table, method=c("range"))
norm_scale1 <- predict(process1, table)
table3 <- melt(norm_scale1, value.name = "Value", id = c("Treatment"), variable.name = "Parameter" )

s1 <- ggplot(table3, aes(x = Parameter, y = Value, fill = factor(Treatment, level = my_order))) +  
  stat_boxplot(geom="errorbar") +
  geom_boxplot() +
  theme_bw() + 
  scale_x_discrete(name = "Treatment") +
  scale_fill_manual(values=c("grey77", "darkorange2", "skyblue3"), name = "Treatment")
s1

ggsave(filename = "WinRhizor_Parameters_SideBySide_Boxplot.pdf", plot = last_plot(), dpi = 600, units = "cm", width = 70, height = 60, scale = 0.5)


#Distribution barplot####
#not normalized
table <- read.csv("DATA_WinRhizor_(distribution).csv", sep=";",
                  header=T)
table2 <- melt(table, value.name = "Value", id = c("Treatment","Replicate"), variable.name = "Diameter_Class" )


Summary_table <- ddply(table2, c("Treatment", "Diameter_Class"), summarise,
                       N    = sum(!is.na(Value)),
                       mean = mean(Value, na.rm=TRUE),
                       sd   = sd(Value, na.rm=TRUE),
                       se   = sd / sqrt(N))

f1 <- ggplot(Summary_table, aes(x = Diameter_Class, y = mean, fill = Treatment)) +  
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9)) +
  scale_x_discrete(name = "Diameter_Class") +
  ggtitle("Distribution barplot") +
  scale_fill_manual(values=c("grey77", "darkorange2", "skyblue3"), name = "Treatment")
f1
ggsave(filename = "WinRhizor_Surface_Distribution_barplot.pdf", plot = last_plot(), dpi = 600, units = "cm", width = 70, height = 60, scale = 0.5)

#normalized per treatment
table3 <- dcast(table2, Diameter_Class + Replicate ~ Treatment, value.var = "Value")

process1 <- preProcess(table3, method=c("range"))
norm_scale1 <- predict(process1, table3)

norm_scale2 <- norm_scale1[,-2]
table4 <- melt(norm_scale2, value.name = "Value", id = c("Diameter_Class"), variable.name = "Treatment" )
table4 <- table4 %>% drop_na(Value)
Summary_table_2 <- ddply(table4, c("Treatment", "Diameter_Class"), summarise,
                       N    = sum(!is.na(Value)),
                       mean = mean(Value, na.rm=TRUE),
                       sd   = sd(Value, na.rm=TRUE),
                       se   = sd / sqrt(N))

f1 <- ggplot(Summary_table_2, aes(x = Diameter_Class, y = mean, fill = Treatment)) +  
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9)) +
  scale_x_discrete(name = "Diameter_Class") +
  ggtitle("Normalized Distribution barplot") +
  scale_fill_manual(values=c("grey77", "darkorange2", "skyblue3"), name = "Treatment")
f1
ggsave(filename = "WinRhizor_Surface_Distribution_Norm_barplot.pdf", plot = last_plot(), dpi = 600, units = "cm", width = 70, height = 60, scale = 0.5)

#Surface/FW boxplot + stat####
table <- read.csv("DATA_WinRhizor_Surface-FW.csv", sep=";",
                  header=T)
name <- "Surface/FW (cm2/g)"

##boxplots & facet wrap ####
my_order<- c("C", "-Fe", "-P")
f1 <- ggplot(table, aes(x= Treatment, Value, fill=Treatment))+  
  stat_boxplot(geom="errorbar", width=0.2)+
  geom_boxplot(width=0.5)+ 
  theme(legend.position = "NONE") + 
  theme_bw() +
  xlab("Treatments") +
  ylab("Surface/FW (cm2/g)") +
  scale_x_discrete(limits=my_order) +
  scale_fill_manual(values=c("darkorange2","skyblue3", "grey77"))
f1
ggsave(filename = "WinRhizor_Surface-FW (cm2-g).pdf", plot = last_plot(), dpi = 600, units = "cm", width = 30, height = 25, scale = 0.5)

##statistics ####

#Assumptions 
## 1. Homogeneity of variances
##Treatment*Time
L_test <- levene_test(table, Value ~ Treatment)
View(L_test)
write.table(L_test, file = "WinRhizor_Surface-FW_Levene_test_results.csv", quote = FALSE, sep = ";")


##2. Normality
##Shapiro-Wilk test for all single treatments
SW_test <- table %>%
  group_by(Treatment) %>%
  shapiro_test(Value)
View(SW_test)
write.table(SW_test, file = "WinRhizor_Surface-FW_ShapiroWilk_test_results.csv", quote = FALSE, sep = ";")

##3. Indipendency
Data are indepent by experimental design!
  

##anova for tukey
OneWay_Anova_Boxplot <- aov(Value ~ Treatment, data = table)

##anova for print
OneWay_Anova_Boxplot2 <- anova(lm(Value ~ Treatment, data = table))

sink("WinRhizor_Surface-FW_OneWay_Anova_Boxplot.csv")
OneWay_Anova_Boxplot2
sink(NULL)

#Tukey
##HSD complete
HSD_Boxplot <- HSD.test(OneWay_Anova_Boxplot, "Treatment")

##HSD groups only
HSD_Boxplot_groups <- as.data.frame(HSD_Boxplot[["groups"]])

sink("WinRhizor_Surface-FW_HSD_Boxplot.csv")
HSD_Boxplot_groups
sink(NULL)
