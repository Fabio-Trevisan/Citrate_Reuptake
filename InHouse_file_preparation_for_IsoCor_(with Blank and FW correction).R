library(dplyr)
library(scales)
library(tidyr)
library(stringr)
library(reshape2)


#AA
#Read CSV ####
df <- read.csv("20230714_AccQ-Tag_AA_peak areas_inhouse_script 2.0.csv", sep=";", header=T)
Class <- "AA"



#OA
#Read CSV ####
df <- read.csv("20230711 3-NPH acids from in-house script.csv", sep=";", header=T)
Class <- "OA"



df <- df[-1,]
df <- data.frame(t(df)) #reverse dataframe 
names(df) <- t(df[1,]) #rename columns
df <- df[-1,]

df[,] <- sapply(df[,], as.numeric) #columns from character to numeric
df <- df[!str_detect(names(df), "QC")]



# Sample/Blank ratio filtering ####
Blank_filtering <- df

y = 3 #blank to sample ratio for filtering
Blank_filtering$Blank <- pmax(Blank_filtering$`010_blank_50% MeOH.cdf`, 
                              Blank_filtering$`011_extraction blank.cdf`) * y #create column with LOD

vector_samples <- names(df) #vector with sample names

Blank_filtering <- as.data.frame(lapply(vector_samples, function(m){
  ifelse(Blank_filtering$Blank > Blank_filtering[,m], 0, Blank_filtering[,m])
})) #create df with value filtered with LOD
names(Blank_filtering) <- vector_samples
row.names(Blank_filtering) <- row.names(df) #rename cols and rows of new df



#clening ratios and means from dataset
Blank_filtering <- Blank_filtering[!str_detect(names(Blank_filtering), "lank")]
df1 <- Blank_filtering[!str_detect(names(Blank_filtering), "std")]



# Sample weight normalization #### 
FW <- read.csv(paste("FreshWeight_", Class, ".csv", sep=""), sep=";", header=T) #get FW or root weight
row.names(FW) <-  unlist(FW[,1]) #rename rows
FW <- as.data.frame(t(FW)) 
FW <- FW[-1,]#remove headings columns/rows
FW[,] <- sapply(FW[,], as.numeric) #set numeric variables
FW <- FW[rep(1,nrow(df1)),] #create df of equal size

FW_norm <- df1/FW #dived the 2 df



# Removal of outlayers####
#FW_norm <- FW_norm[,!(names(FW_norm) %in% c("16_M2","48_P6"))]



#Dataset cleaning and save ####
FW_norm$metabolite <- row.names(FW_norm)
row.names(FW_norm) <- NULL
FW_norm$derivative <- 0
FW_norm$isotopologue <- 0
FW_norm$resolution <- 0 #add rown needed by IsoCor

df2 <- melt(FW_norm, value.name = "area", variable.name = "sample", 
            id = c("metabolite", "derivative", "isotopologue", "resolution")) #melt columns

IsoCor <- read.csv(paste("IsoCor_FileForm_", Class, ".csv", sep=""), sep=";", header=T) #get FW or root weight
df2$metabolite <- IsoCor$metabolite
row.names(df2) <- NULL
df2$derivative <- IsoCor$derivative
df2$isotopologue <- IsoCor$isotopologue
df2$resolution <- IsoCor$resolution #add rown needed by IsoCor

if (Class == "OA") {write.table(df2, file = "20230711 3-NPH acids from in-house script_IsoCor.csv", sep =";", row.names=FALSE)} else {write.table(df2, file = "20230714_AccQ-Tag_AA_peak areas_inhouse_script 2.0_IsoCor.csv", sep =";", row.names=FALSE)}

if (Class == "OA") {write.table(df2, file = "20230711 3-NPH acids from in-house script_IsoCor.tsv", sep ="\t", quote=FALSE, row.names=FALSE)} else {write.table(df2, file = "20230714_AccQ-Tag_AA_peak areas_inhouse_script 2.0_IsoCor.tsv", sep ="\t", quote=FALSE, row.names=FALSE)}
