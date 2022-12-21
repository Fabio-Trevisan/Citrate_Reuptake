library(plyr) 
library(dplyr)
library(reshape2)
library(caret)
library(ggplot2)


#Read csv ####
table <- read.csv("DATA_Citrate_ICP.csv", sep=";",
                  header=T)


#Re-arrange table and mean value calculation ####
table2 <- melt(data = table, 
               id.vars = c("Tissue", "Treatment", "Concentration"), 
               variable.name = "Element", 
               value.name = "ppm"
)

table3 <- dcast(table2, 
                Tissue + Treatment + Concentration ~Element, 
                mean,
                value.var = "ppm") 


#Normalization and centration for each element####
#between 0 and 1 (according to range)
process1 <- preProcess(table3, method=c("range"))
norm_scale1 <- predict(process1, table3)

norm_scale <- melt(norm_scale1, id = c("Tissue","Treatment", "Concentration"))

a <- ggplot(norm_scale, aes(Treatment, variable, size=value, colour = Treatment)) +
  geom_point() + 
  facet_wrap(~Tissue + Concentration, ncol = 2)
a

b <- ggplot(norm_scale, aes(Treatment, variable, size=value, colour = Treatment)) +
  geom_point() + 
  facet_wrap(~Tissue, ncol = 2)
b

ggsave(filename = "ICP_Citrate_Bubbleplot.pdf", 
       plot = last_plot(), 
       dpi = 600, 
       units = "cm", 
       width = 70, 
       height = 70, 
       scale = 0.2)

c <- ggplot(norm_scale, aes(Concentration, variable, size=value, colour = Concentration)) +
  geom_point() + 
  facet_wrap(~Tissue + Treatment, ncol = 2)
c

ggsave(filename = "ICP_Citrate_Bubbleplot_Conc.Comparison.pdf", 
       plot = last_plot(), 
       dpi = 600, 
       units = "cm", 
       width = 45, 
       height = 80, 
       scale = 0.35)


