library(dplyr)
library(scales)
library(tidyr)
library(stringr)
library(reshape2)

library(tidyverse)
library(ggpubr)


# AA ####
#Read CSV
df <- read.csv("20230714_AccQ-Tag_AA_peak areas_inhouse_script 2.0.csv", sep=";", header=T)
IsoCor <- read.csv("20230714_AccQ-Tag_AA_peak areas_inhouse_script 2.0_IsoCor_res.tsv", sep="\t", header=T)
Class <- "AA"
Class_2 <- "A"
Std.Conc <- c(1,10,100)
m0 <- c("Ala_260", 
        "Arg_345", 
        "Asn_303", 
        "Asp_304",
        "Citrulline_346",
        "Gln_317", 
        "Gly_246", 
        "His_326", 
        "Ile_302", 
        "Leu_302", 
        "Lys_487", 
        "Met_320", 
        "Phe_336",
        "Pro_286",
        "Ser_276",
        "Thr_290",
        "Trp_375", 
        "Tyr_352", 
        "Val_288",
        "GABA_274",
        "Glu_318",
        "Ornithine_473") #list of m+0 isotopologues for filtering

m1 <- c("Alanine", 
        "Arginine", 
        "Asparagine", 
        "Aspartate",
        "Citrulline",
        "Glutamine", 
        "Glycine", 
        "Histidine", 
        "Isoleucine", 
        "Leucine", 
        "Lysine", 
        "Methionine", 
        "Phenylalanine",
        "Proline",
        "Serine",
        "Threonine",
        "Tryptophan", 
        "Tyrosine", 
        "Valine",
        "GABA",
        "Glutamate",
        "Ornithine")

# OA ####
#Read CSV
df <- read.csv("20230711 3-NPH acids from in-house script.csv", sep=";", header=T)
IsoCor <- read.csv("20230711 3-NPH acids from in-house script_IsoCor_res.tsv", sep="\t", header=T)
Class <- "OA"
Class_2 <- "O"
Std.Conc <- c(1,10,50,100)
m0 <- c("lactate_224", 
        "pyruvate_357", 
        "fumarate_385", 
        "succinate_387", 
        "malate_403", 
        "isocitrate_443", 
        "citrate_443", 
        "oxaloacetate_536", 
        "X2.oxoglutarate_550", 
        "cis.aconitate_578", 
        "glycerate_240", 
        "glycolate_210",
        "X2.hydroxyglutarate_417",
        "oxalate_359",
        "itaconate_399") #list of m+0 isotopologues for filtering

m1 <- c("Lactate",
        "Pyruvate",
        "Fumarate",
        "Succinate",
        "Malate",
        "Isocitrate",
        "Citrate",
        "Oxaloacetate",
        "Oxoglutarate",
        "Aconitate",
        "Glycerate",
        "Glycolate",
        "Hydroxyglutarate",
        "Oxalate",
        "Itaconate")


# df curation ####
df <- df[-1,]
df <- data.frame(t(df)) #reverse dataframe 
names(df) <- t(df[1,]) #rename columns
df <- df[-1,]

df[,] <- sapply(df[,], as.numeric) #columns from character to numeric
df <- df[!str_detect(names(df), "QC")] #remove QC
df$Delete <- row.names(df)

df <- df[df$Delete %in% m0, ]
df <- df[!str_detect(names(df), "Delete")] #remove QC


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
Blank_filtered <- Blank_filtering[!str_detect(names(Blank_filtering), "lank")] #cleaning blanks from df



# Concentration calculation ####
std <- Blank_filtered[!str_detect(names(Blank_filtered), Class_2)] #select only std
vector_metabolite <- row.names(std) #create vector with metabolites names
t_std <- as.data.frame(t(std)) #transpose df
t_std$Std_conc <- Std.Conc #create columns with standard concentration uM
n <- ncol(t_std)
Regression <- lapply(t_std[1:n], function(m){
  lm(m ~ Std_conc, data = t_std)
}) #create linear regression model

Coefficients <- lapply(vector_metabolite, function(m){
  as.data.frame(Regression[[m]][["coefficients"]])
}) #extract  intercept and slope
names(Coefficients) <- vector_metabolite #add names
for(i in vector_metabolite) {
      names(Coefficients[[i]]) <- "Coefficients"
} #rename column to Coefficients


## Transform list in dataframe ####
Coefficients_2 <- lapply(vector_metabolite, function(m){
  as.data.frame(t(Coefficients[[m]][["Coefficients"]]))
}) #transpose coefficients
names(Coefficients_2) <- vector_metabolite #add names
for(i in vector_metabolite) {
  names(Coefficients_2[[i]]) <- c("Intercept", "Slope")
} #rename columns to Intercept and Slope respectively

Coefficients_merged <- do.call(rbind.data.frame, Coefficients_2) #merged coefficients into one df


## Conversion area to concentration ####
x <- ncol(Blank_filtered)-1 #counter adapting df size
Intercept <- as.data.frame(Coefficients_merged$Intercept) #create df only with intercept values
Intercept <- cbind(Intercept, rep(Intercept[1],x)) #adapt df size
Slope <- as.data.frame(Coefficients_merged$Slope) #create df only with slope values
Slope <- cbind(Slope, rep(Slope[1],x)) #adapt df size

Concentration <- (Blank_filtered - Intercept) / Slope #calculate concentration according to x=(y-q)/m inverse of (y=mx+q)
#save concentration table
Concentration <- Concentration[!str_detect(names(Concentration), "std")] #remove stds
Concentration$metabolite <- row.names(Concentration) #create variable metabolites
row.names(Concentration) <- NULL #remove row names
Concentration$metabolite <- m1 #rename properly

Conversion <- melt(Concentration, variable.name = "sample", id.vars = "metabolite", value.name = "concentration")

# Import isocor results ####
## cleaning
IsoCor <- IsoCor[,-c(3,5:7,9:10)] #remove un-useful columns
IsoCor <- IsoCor[!str_detect(IsoCor$sample, "blank"),] #remove blank samples
IsoCor <- IsoCor[!str_detect(IsoCor$sample, "std"),] #remove standards
IsoCor <- IsoCor[!str_detect(IsoCor$sample, "QC"),] #remove QCs

IsoCor <- IsoCor[!str_detect(IsoCor$isotopologue, "0"),] #remove m+0
IsoCor <- IsoCor[!str_detect(IsoCor$metabolite, "Norvaline"),] #remove m+0

IsoCor$isotopologue <- as.factor(IsoCor$isotopologue)
Isotopologues <- levels(IsoCor$isotopologue)

## Add concentration to IsoCor df
Conversion <- Conversion
IsoCor <- IsoCor
Conversion$sample <- as.factor(Conversion$sample)
Conversion$metabolite <- as.factor(Conversion$metabolite)
IsoCor$sample <- as.factor(IsoCor$sample)
IsoCor$metabolite <- as.factor(IsoCor$metabolite)

IsoCor$concentration <- 0

for(i in 1:nrow(Conversion)) {
  IsoCor$concentration[IsoCor$sample == Conversion[i,]$sample & IsoCor$metabolite == Conversion[i,]$metabolite] <- Conversion[i,]$concentration
}


## Calculate 13C concentration for each isotopologue (Excess AA+n * [AA] * n = [13C AA+n])




lv <- union(levels(two$sample), levels(two$metabolite))

one$sample <- factor(one$sample, levels = lv)
one$metabolite <- factor(one$metabolite, levels = lv)
two$sample <- factor(two$sample, levels = lv)
two$metabolite <- factor(two$metabolite, levels = lv)



##Calculate concentration of isotopologues with the following formula: 
###Excess AA+n * [AA] = [AA+n]
###or
###Excess OA+n * [OA] = [OA+n]
### Where  Excess OA+n = %13C Enrichment 










##Sum concentration of isotopologues * their mass compared with the monoisotopic form (e.g. m+2 will be multiblied by 2)
##and molecules (allla AA and all OA)

##Account for dilution 
### AA 1mL extraction down to 100microL
### OA 1mL extraction down to 50microL 
### FW already used before IsoCor (but only %13C enrichment was used from IsoCor, since the concentration was calculated on the raw data FW normalisation again)
### Divide total concentration by FW 

### [C]microM * Volume(0.1 or 0.05mL) / FW (mg) * 10 (FW/DW) = micromol 13C / g DW
### check FW / DW ratio
### Combine AA and OA



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











# IRMS data ####
#Read CSV ####
df <- read.csv("---.csv", sep=";", header=T)

#convert stable isotope ratio in carbon content and 13C carbon concentration 
#via mass balance / mm of citrate?
#or
#via stable carbon isotope ratio knowing the concentration of C in the sample? 