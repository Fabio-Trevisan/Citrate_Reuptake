library(plyr)
library(dplyr)
library(scales)
library(tidyr)
library(stringr)
library(reshape2)

library(tidyverse)
library(ggpubr)

# CSIA (run first script for AA and then for OA)
## OA ####
#Read CSV
df <- read.csv("20230711 3-NPH acids from in-house script.csv", sep=";", header=T)
IsoCor <- read.csv("20230711 3-NPH acids from in-house script_IsoCor_res.tsv", sep="\t", header=T)
Class <- "OA"
Class_2 <- "O"
Dilution <- 0.1 #mL  or 100uM used for resuspension after extraction
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


## AA ####
#Read CSV
df <- read.csv("20230714_AccQ-Tag_AA_peak areas_inhouse_script 2.0.csv", sep=";", header=T)
IsoCor <- read.csv("20230714_AccQ-Tag_AA_peak areas_inhouse_script 2.0_IsoCor_res.tsv", sep="\t", header=T)
Class <- "AA"
Class_2 <- "A"
Dilution <- 0.05 #mL  or 50uM used for resuspension after extraction
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

## df curation ####
df <- df[-1,]
df <- data.frame(t(df)) #reverse dataframe 
names(df) <- t(df[1,]) #rename columns
df <- df[-1,]

df[,] <- sapply(df[,], as.numeric) #columns from character to numeric
df <- df[!str_detect(names(df), "QC")] #remove QC
df$Delete <- row.names(df)

df <- df[df$Delete %in% m0, ]
df <- df[!str_detect(names(df), "Delete")] #remove QC


## Sample/Blank ratio filtering ####
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



## Concentration calculation ####
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


### Transform list in dataframe ####
Coefficients_2 <- lapply(vector_metabolite, function(m){
  as.data.frame(t(Coefficients[[m]][["Coefficients"]]))
}) #transpose coefficients
names(Coefficients_2) <- vector_metabolite #add names
for(i in vector_metabolite) {
  names(Coefficients_2[[i]]) <- c("Intercept", "Slope")
} #rename columns to Intercept and Slope respectively

Coefficients_merged <- do.call(rbind.data.frame, Coefficients_2) #merged coefficients into one df


### Conversion area to concentration ####
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

### Import isocor results ####
#### cleaning ####
IsoCor <- IsoCor[,-c(3,5:7,9:10)] #remove un-useful columns
IsoCor <- IsoCor[!str_detect(IsoCor$sample, "blank"),] #remove blank samples
IsoCor <- IsoCor[!str_detect(IsoCor$sample, "std"),] #remove standards
IsoCor <- IsoCor[!str_detect(IsoCor$sample, "QC"),] #remove QCs

IsoCor <- IsoCor[!str_detect(IsoCor$isotopologue, "0"),] #remove m+0
IsoCor <- IsoCor[!str_detect(IsoCor$metabolite, "Norvaline"),] #remove m+0

IsoCor$sample <- gsub("^.{0,4}", "", IsoCor$sample)
IsoCor$sample <- gsub(".{4}$", "", IsoCor$sample)
Conversion$sample <- gsub("^.{0,4}", "", Conversion$sample)
Conversion$sample <- gsub(".{4}$", "", Conversion$sample)

#### Add concentration to IsoCor df ####
Conversion$metabolite <- as.factor(Conversion$metabolite) #convert character or numeric variable into factor
IsoCor$metabolite <- as.factor(IsoCor$metabolite)

IsoCor$concentration <- 0 #create variable in IsoCor df

for(i in 1:nrow(Conversion)) {
  IsoCor$concentration[IsoCor$sample == Conversion[i,]$sample & IsoCor$metabolite == Conversion[i,]$metabolite] <- Conversion[i,]$concentration
} #add concentration con conversion df in IsoCor df


## Calculate 13C concentration ####
IsoCor$Carbon13Conc <- 
  IsoCor$isotopologue * IsoCor$isotopologue_fraction * IsoCor$concentration # Calculate 13C concentration (Excess AA+n * [AA] * n = [13C AA+n]) 

IsoCor[is.na(IsoCor)] <- 0 #substitute NA with 0
IsoCor$Carbon13Conc <- ifelse(IsoCor$Carbon13Conc < 0, 0, IsoCor$Carbon13Conc)#remove negative values

IsoCor$sample_number <- substr(IsoCor$sample,2,4) #select sample number
IsoCor$sample_number <- as.numeric(IsoCor$sample_number) #convert character to number

C13Content <- dcast(IsoCor, metabolite ~ sample_number, sum,  value.var = "Carbon13Conc") #cast dataset (long to wide and sum [C] of isotopologues)
row.names(C13Content) <- C13Content$metabolite #rename rows df
C13Content <- C13Content[,-1] #remove first row

### Import dataset FW ####
FW <- read.csv(paste("FreshWeight_", Class, ".csv", sep=""), sep=";", header=T) #get FW (or root weight) weighted for LC-MS analysis approx. 15mg
row.names(FW) <-  unlist(FW[,1]) #rename rows
FW <- as.data.frame(t(FW))#transpose df
FW <- FW[-1,]#remove headings columns/rows
FW[,] <- sapply(FW[,], as.numeric) #set numeric variables
FW <- FW[rep(1,nrow(C13Content)),] #create df of equal size


### Calculate 13C content ####
C13Content2 <- C13Content * Dilution / FW * 20 #[C]uM * Volume(0.1 or 0.05mL) / FW (mg) * 20 (FW/DW) = umol 13C / g DW

if (Class == "OA") C13ContentOA <- as.data.frame(colSums(C13Content2)) else C13ContentAA <- as.data.frame(colSums(C13Content2))

## REPEAT TILL HERE WITH BOTH DATASET "AA" AND "OA" ####
## Combine AA and OA ####
C13ContentAA <- as.data.frame(C13ContentAA[-39,]) #remove samples where only AA was present & keep only common once
C13ContentAA <- as.data.frame(C13ContentAA[-22,])

CSIA_13C <- C13ContentAA + C13ContentOA #sum AA & OA df


## Import Time-Treatment-Labeling factors file
Treatment_factors <- read.csv(paste("Treatment_factors_", "OA", ".csv", sep=""), sep=";", header=T) #load file with treatment_factors
Treatment_factors <- Treatment_factors[1:nrow(CSIA_13C),]

Treatment_factors$C13 <- CSIA_13C$`C13ContentAA[-22, ]`

## Uniformation and vectors creation
Treatment_factors$Treatment <- factor(Treatment_factors$Treatment)
Treatment_factors$Time <- factor(Treatment_factors$Time)
Treatment_factors$Labeling <- factor(Treatment_factors$Labeling) #convert character and numebers to factors

Treatment_factors_L <- Treatment_factors %>% filter(str_detect(Labeling, "L"))
Treatment_factors_L <- Treatment_factors_L[,-4]
Treatment_factors_L <- Treatment_factors_L[,-1]
colnames(Treatment_factors_L)[3] ="Value"

### Summary table ####
Summary_table_CSIA <- ddply(Treatment_factors_L, c("Time", "Treatment"), summarise,
                       N    = sum(!is.na(Value)),
                       mean = mean(Value, na.rm=TRUE),
                       sd   = sd(Value, na.rm=TRUE),
                       se   = sd / sqrt(N))





# BSIA ####
BSIA <- read.csv("DATA_STD_MassBalance_Citrate_RECOVERY.csv", sep=";",
                  header=T)

BSIA$Species_Tissue <- as.factor(BSIA$Species_Tissue)
BSIA$Time <- as.factor(BSIA$Time)
BSIA$Treatment <- as.factor(BSIA$Treatment)

BSIA$Value <- BSIA$Value/13*1000
BSIA <- BSIA[,-3]

## Summary table
Summary_table_BSIA <- ddply(BSIA, c("Time", "Treatment"), summarise,
                       N    = sum(!is.na(Value)),
                       mean = mean(Value, na.rm=TRUE),
                       sd   = sd(Value, na.rm=TRUE),
                       se   = sd / sqrt(N))


# Citrate 2+ and tot-isotopologues ####
Citrate <- IsoCor[IsoCor$metabolite =="Citrate",]
Citrate_M2 <- IsoCor[IsoCor$metabolite =="Citrate" & IsoCor$isotopologue =="2",]

Citrate <- dcast(Citrate, metabolite ~ sample_number, sum,  value.var = "Carbon13Conc") #cast dataset (long to wide and sum [C] of isotopologues)
row.names(Citrate) <- Citrate$metabolite #rename rows df
Citrate <- Citrate[,-1] #remove first row

Citrate_M2 <- dcast(Citrate_M2, metabolite ~ sample_number, sum,  value.var = "Carbon13Conc") #cast dataset (long to wide and sum [C] of isotopologues)
row.names(Citrate_M2) <- Citrate_M2$metabolite #rename rows df
Citrate_M2 <- Citrate_M2[,-1] #remove first row

## Import dataset FW ####
FW <- read.csv(paste("FreshWeight_", Class, ".csv", sep=""), sep=";", header=T) #get FW (or root weight) weighted for LC-MS analysis approx. 15mg
row.names(FW) <-  unlist(FW[,1]) #rename rows
FW <- as.data.frame(t(FW))#transpose df
FW <- FW[-1,]#remove headings columns/rows
FW[,] <- sapply(FW[,], as.numeric) #set numeric variables

## Calculate 13C content ####
Citrate <- Citrate * Dilution / FW * 20 #[C]uM * Volume(0.1 or 0.05mL) / FW (mg) * 20 (FW/DW) = umol 13C / g DW
Citrate_M2 <- Citrate_M2 * Dilution / FW * 20 #[C]uM * Volume(0.1 or 0.05mL) / FW (mg) * 20 (FW/DW) = umol 13C / g DW

Citrate <- as.data.frame(t(Citrate))
Citrate_M2 <- as.data.frame(t(Citrate_M2))

## Import Time-Treatment-Labeling factors file
Treatment_factors <- read.csv(paste("Treatment_factors_", "OA", ".csv", sep=""), sep=";", header=T) #load file with treatment_factors
Citrate_TF <- Treatment_factors[1:117,]
Citrate_M2_TF <- Treatment_factors[1:117,]

Citrate_TF$C13 <- Citrate$Citrate
Citrate_M2_TF$C13 <- Citrate_M2$Citrate


## Uniformation and vectors creation
Citrate_TF$Treatment <- factor(Citrate_TF$Treatment)
Citrate_TF$Time <- factor(Citrate_TF$Time)
Citrate_TF$Labeling <- factor(Citrate_TF$Labeling) #convert character and numebers to factors

Citrate_M2_TF$Treatment <- factor(Citrate_M2_TF$Treatment)
Citrate_M2_TF$Time <- factor(Citrate_M2_TF$Time)
Citrate_M2_TF$Labeling <- factor(Citrate_M2_TF$Labeling) #convert character and numebers to factors

Citrate_final <- Citrate_TF %>% filter(str_detect(Labeling, "L"))
Citrate_final <- Citrate_final[,-4]
Citrate_final <- Citrate_final[,-1]
colnames(Citrate_final)[3] ="Value"

Citrate_M2_final <- Citrate_M2_TF %>% filter(str_detect(Labeling, "L"))
Citrate_M2_final <- Citrate_M2_final[,-4]
Citrate_M2_final <- Citrate_M2_final[,-1]
colnames(Citrate_M2_final)[3] ="Value"

## Summary table ####
Summary_table_Citrate <- ddply(Citrate_final, c("Time", "Treatment"), summarise,
                            N    = sum(!is.na(Value)),
                            mean = mean(Value, na.rm=TRUE),
                            sd   = sd(Value, na.rm=TRUE),
                            se   = sd / sqrt(N))

Summary_table_Citrate_M2 <- ddply(Citrate_M2_final, c("Time", "Treatment"), summarise,
                               N    = sum(!is.na(Value)),
                               mean = mean(Value, na.rm=TRUE),
                               sd   = sd(Value, na.rm=TRUE),
                               se   = sd / sqrt(N))





# Combine CSIA & BSIA & Citrate ####
Summary_table_combined <- rbind(Summary_table_BSIA, Summary_table_CSIA, Summary_table_Citrate, Summary_table_Citrate_M2)
Table_combined <- rbind(BSIA, Treatment_factors_L, Citrate_final, Citrate_M2_final)


## Save table (umol 13C / g DW)  [assuming FW = 5% DW according to the data we have] ####
write.table(Summary_table_combined, file = "CSIA_Recovery_summary.csv", sep =";", row.names=FALSE)

write.table(Table_combined, file = "CSIA_Recovery_raw.csv", sep =";", row.names=FALSE)
