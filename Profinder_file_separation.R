
df <- read.csv("20230525_export_AccQ-Tag_AA.csv", sep = ";" , head = F)

df2 <- split(df, df$V1)

df2 <- df2[-1]

lapply(names(df2), function(x) {
  x1 <- df2[[x]]
  write.table(x1, file = paste(x, ".csv", sep =""))
})
