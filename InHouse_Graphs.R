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
library(ggbreak)


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
df["mean_enrichment"][df["mean_enrichment"]==0] <- NA #raplace meanenrichment 0 values with 1/10 ot 1/100 of min values
x <- min(df$mean_enrichment, na.rm = T)/10 #calculate min values/10
df[is.na(df)] <-  runif(sum(is.na(df)), min = x-(x/10), max = x+(x/10))



# Vector Significant Molecules ####
OA_Significant_Molecules <- "Aconitate|Citrate|Hydroxyglutarate|Isocitrate|Oxoglutarate"
#Pyruvate, Succinare, Oxalate and Oxaloacetate P<0.05 but removed because non-sense 

AA_Significant_Molecules <- "Aspartate|Citrulline|GABA|Glutamate|Glutamine"
#Lysine and Alanine 0.05<P<0.10 



# Summary table ####
Summary_table <- ddply(df, c("metabolite", "Time", "Labeling", "Treatment"), summarise,
                       N    = sum(!is.na(mean_enrichment)),
                       mean = mean(mean_enrichment, na.rm=TRUE),
                       sd   = sd(mean_enrichment, na.rm=TRUE),
                       se   = sd / sqrt(N))

Summary_table_L <- Summary_table %>% filter(str_detect(Labeling, "L"))



# Barplot ####
##   Significant only 
if (Class == "AA") {
  Summary_table_L_Significant <- Summary_table_L %>% 
    filter(str_detect(metabolite, gsub('["]', '', AA_Significant_Molecules)))
} else { 
  Summary_table_L_Significant <- Summary_table_L %>% 
    filter(str_detect(metabolite, gsub('["]', '', OA_Significant_Molecules)))
}

f3 <- ggplot(Summary_table_L_Significant, aes(x = Time, y = mean, fill = Treatment)) + 
  geom_bar(stat = "identity", position = 'dodge', width = 0.8) + 
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, group = Treatment), width = 0.5 , position = position_dodge(.8)) +
  scale_fill_manual(values=c("grey77", "darkorange2", "skyblue3")) +
  theme_bw() + 
  scale_y_continuous(labels = percent)
f4 <- f3 + facet_wrap(~metabolite, scales="fixed") + 
  ylab("% Enrichment") + 
  xlab("Time (min)") 
f4

ggsave(filename = paste("InHouse_", Class, "_Barplot_Significant.pdf", sep=""), plot = f4, dpi = 600, units = "cm", width = 80, height = 60, scale = 0.35)



# Scatterplot ####
##   Significant only 
Summary_table_L$Time <- as.character(Summary_table_L$Time)
Summary_table_L$Time <- as.numeric(Summary_table_L$Time)
if (Class == "AA") {
  Summary_table_L_Significant <- Summary_table_L %>% 
    filter(str_detect(metabolite, gsub('["]', '', AA_Significant_Molecules)))
} else { 
  Summary_table_L_Significant <- Summary_table_L %>% 
    filter(str_detect(metabolite, gsub('["]', '', OA_Significant_Molecules)))
}

f3 <- ggplot(Summary_table_L_Significant, aes(x = Time, y = mean, group = Treatment, colour = Treatment)) + 
  geom_line(aes(group = Treatment)) + 
  geom_point(aes(shape = Treatment)) + 
  scale_shape_manual(values = c(15:18)) +
  scale_color_manual(values=c("grey77", "darkorange2", "skyblue3"))+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, group = Treatment), width = 5) +
  theme_bw() + 
  scale_y_continuous(labels = percent, breaks = seq(0, 0.07, 0.005)) +
  scale_x_continuous(breaks=seq(0,120,15))
f4 <- f3 + facet_wrap(~metabolite, scales="free") + 
  ylab("% Enrichment") + 
  xlab("Time (min)") 
f4

ggsave(filename = paste("InHouse_", Class, "_Scatterplot_Significant_2.pdf", sep=""), plot = f4, dpi = 600, units = "cm", width = 80, height = 60, scale = 0.35)



# SideBySide Boxlot + trendline ####
df2 <- df %>% filter(str_detect(Labeling, "L"))
if (Class == "AA") {
  df3 <- df2 %>% 
    filter(str_detect(metabolite, gsub('["]', '', AA_Significant_Molecules)))
} else { 
  df3 <- df2 %>% 
    filter(str_detect(metabolite, gsub('["]', '', OA_Significant_Molecules)))
}

trendline <- geom_smooth(aes(group = Treatment, color = Treatment, fill = Treatment), method=lm, alpha = 0.1, linetype="dashed") 

f5 <- ggplot(df3, aes(x = factor(Time), y = mean_enrichment, fill = Treatment)) +  
  trendline +
  stat_boxplot(geom="errorbar") +
  geom_boxplot() +
  theme_bw() + 
  scale_fill_manual(values=c("grey77","darkorange2", "skyblue3")) +
  scale_color_manual(values=c("grey77","darkorange2", "skyblue3")) +
  scale_y_continuous(labels = percent) +
  stat_regline_equation()

f6 <- f5 + facet_wrap(~metabolite, scales="fixed")+
  ylab("% Enrichment") + 
  xlab("Time (min)") 
f6

ggsave(filename = paste("InHouse_", Class, "_BoxPlot_Significant.pdf", sep=""), plot = f6, dpi = 600, units = "cm", width = 80, height = 60, scale = 0.35)

# Histogram ISOTOPOLOGUES ####
##   Dataset preparation ####
#cleaning
df <- table[,-c(3,5:7,9:10)] #remove un-useful columns
df <- df[!str_detect(df$sample, "blank"),] #remove blank samples
df <- df[!str_detect(df$sample, "std"),] #remove standards
df <- df[!str_detect(df$sample, "QC"),] #remove QCs
row.names(df) <- NULL

#Time-Treatment-Labeling 
Treatment_factors <- read.csv(paste("Treatment_factors_", Class, "_isotopologues.csv", sep=""), sep=";", header=T) #load file with treatment_factors
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
df["isotopologue_fraction"][df["isotopologue_fraction"]==0] <- NA #raplace meanenrichment 0 values with 1/10 ot 1/100 of min values
x <- min(df$isotopologue_fraction, na.rm = T)/10 #calculate min values/10
df[is.na(df)] <-  runif(sum(is.na(df)), min = x-(x/10), max = x+(x/10))



##   Summary table ####
Summary_table <- ddply(df, c("metabolite", "Time", "Labeling", "Treatment", "isotopologue"), summarise,
                       N    = sum(!is.na(isotopologue_fraction)),
                       mean = mean(isotopologue_fraction, na.rm=TRUE),
                       sd   = sd(isotopologue_fraction, na.rm=TRUE),
                       se   = sd / sqrt(N))

Summary_table_L <- Summary_table %>% filter(str_detect(Labeling, "L"))



##    Graph ####
OA_Significant_Molecules <- "Aconitate|Citrate|Hydroxyglutarate|Isocitrate|Oxoglutarate"
#Pyruvate, Succinare, Oxalate and Oxaloacetate P<0.05 but removed because non-sense 

AA_Significant_Molecules <- "Aspartate|Citrulline|GABA|Glutamate|Glutamine"
#Lysine and Alanine 0.05<P<0.10 

if (Class == "AA") {
  Summary_table_L_Significant <- Summary_table_L %>% 
    filter(str_detect(metabolite, gsub('["]', '', AA_Significant_Molecules)))
} else { 
  Summary_table_L_Significant <- Summary_table_L %>% 
    filter(str_detect(metabolite, gsub('["]', '', OA_Significant_Molecules)))
}
if (Class == "AA") {
  a <- c(0.25, 0.79)
  b <- 1
} else { 
  a <- c(0.15, 0.84)
  b <- 1.3
}

f7 <- ggplot(Summary_table_L_Significant, aes(x = isotopologue, y = mean, fill = Time)) + 
  geom_bar(stat = "identity", position = 'dodge', width = 0.8) + 
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, group = Time), width = 0.5, position = position_dodge(0.8)) +
  scale_fill_manual(values=c("grey77", "skyblue3", "blue", "navy")) +
  theme_bw() + 
  scale_x_continuous(breaks=seq(0,6,1))

f8 <- f7 +  facet_wrap(~metabolite+Treatment, scales="fixed", ncol = 3) +
  ylab("% Enrichment") + 
  xlab("Isotopologue") +
  scale_y_break(a,  scales = b) +
  scale_y_continuous(breaks = seq(0, 1 , 0.05), labels = percent) 
f8

ggsave(filename = paste("InHouse_", Class, "_Barplot_Isotopologues.pdf", sep=""), plot = f8, dpi = 600, units = "cm", width = 100, height = 150, scale = 0.3)
