library(readr)
library(ggpubr)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(plyr)
library(dplyr)



#Scatter line + error bars (STD.Err) ####

table2 <- read.csv("DATA_13C_Citrate.csv", sep=";",
                   header=T)


#table2 <- table2 %>% drop_na(value)


Summary_table <- ddply(table2, c("Treatment", "Time", "Species_Tissue"), summarise,
                       N    = sum(!is.na(Value)),
                       mean = mean(Value, na.rm=TRUE),
                       sd   = sd(Value, na.rm=TRUE),
                       se   = sd / sqrt(N))
Summary_table


f2 <- ggplot(Summary_table, aes(x = Time, y = mean, group = Treatment, colour = Treatment)) + 
  geom_line(aes(group = Treatment)) + 
  geom_point(aes(shape = Treatment)) + 
  scale_shape_manual(values = c(15:18)) +
  scale_color_manual(values=c("grey77", "darkorange2", "skyblue3"))+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, group = Treatment), width = 0.5) +
  theme_bw() + 
  scale_y_continuous(labels = scientific) +
  scale_x_continuous(breaks=seq(0,17,1))
f2 + facet_wrap(~Species_Tissue, scales="free", ncol = 2) + 
  ylab("Delta 13C") + 
  xlab("Time (Days)") 


ggsave(filename = "13C_Scatter-lines_errorbars(bw)_Citrate.pdf", plot = last_plot(), dpi = 600, units = "cm", width = 70, height = 80, scale = 0.5)


write.table(Summary_table, "13C_Citrate_summary_statistics.csv", quote = FALSE, sep = ";")


#######

table2 <- read.csv("DATA_13C_Citrate_2.csv", sep=";",
                   header=T)


#table2 <- table2 %>% drop_na(value)


Summary_table <- ddply(table2, c("Treatment", "Time", "Species", "Tissue"), summarise,
                       N    = sum(!is.na(Value)),
                       mean = mean(Value, na.rm=TRUE),
                       sd   = sd(Value, na.rm=TRUE),
                       se   = sd / sqrt(N))
Summary_table
 
   
f2 <- ggplot(Summary_table, aes(x = Time, y = mean, group = Species, colour = Species)) + 
  geom_line(aes(group = Species)) + 
  geom_point(aes(shape = Species)) + 
  scale_shape_manual(values = c(15:18)) +
  scale_color_manual(values=c("grey77", "darkorange2", "skyblue3"))+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, group = Species), width = 0.5) +
  theme_bw() + 
  scale_y_continuous(labels = scientific) +
  scale_x_continuous(breaks=seq(0,17,1))
f2 + facet_wrap(~Treatment + Tissue, scales="free", ncol = 2) + 
  ylab("Delta 13C") + 
  xlab("Time (Days)") 

ggsave(filename = "13C_Scatter-lines_errorbars_conc-comparison_Citrate.pdf", plot = last_plot(), dpi = 600, units = "cm", width = 70, height = 80, scale = 0.5)

