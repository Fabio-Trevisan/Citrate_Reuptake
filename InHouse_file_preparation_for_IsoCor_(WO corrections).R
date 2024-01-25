library(readr)
library(ggpubr)
library(tidyverse)
library(plyr)
library(dplyr)
library(rstatix)
library(purrr)
library(agricolae)
library(reshape2)

#OA
df <- read.csv("20230711 3-NPH acids from in-house script.csv", sep=";",
               header=T)
Class = "OA"

#AA
df <- read.csv("20230714_AccQ-Tag_AA_peak areas_inhouse_script 2.0.csv", sep=";",
               header=T)
Class = "AA"

#Script
df <- df[-1,]
df <- data.frame(t(df)) #reverse dataframe 
names(df) <- t(df[1,]) #rename columns
df <- df[-1,]

df <- df[!str_detect(names(df), "lank")]
df <- df[!str_detect(names(df), "std")]
df <- df[!str_detect(names(df), "QC")]

df$metabolite <- row.names(df)
row.names(df) <- NULL
df$derivative <- 0
df$isotopologue <- 0
df$resolution <- 0 #add rown needed by IsoCor

df1 <- melt(df, value.name = "area", variable.name = "sample", 
            id = c("metabolite", "derivative", "isotopologue", "resolution")) #melt columns

IsoCor <- read.csv(paste("IsoCor_FileForm_", Class, ".csv", sep=""), sep=";", header=T) #get FW or root weight
df1$metabolite <- IsoCor$metabolite
row.names(df1) <- NULL
df1$derivative <- IsoCor$derivative
df1$isotopologue <- IsoCor$isotopologue
df1$resolution <- IsoCor$resolution #add rown needed by IsoCor


if (Class == "OA") {write.table(df1, file = "20230711 3-NPH acids from in-house script_IsoCor.csv", sep =";", row.names=FALSE)} else {write.table(df1, file = "20230714_AccQ-Tag_AA_peak areas_inhouse_script 2.0_IsoCor.csv", sep =";", row.names=FALSE)}

if (Class == "OA") {write.table(df1, file = "20230711 3-NPH acids from in-house script_IsoCor.tsv", sep ="\t", quote=FALSE, row.names=FALSE)} else {write.table(df1, file = "20230714_AccQ-Tag_AA_peak areas_inhouse_script 2.0_IsoCor.tsv", sep ="\t", quote=FALSE, row.names=FALSE)}
