library(tcltk)
library(tkrplot)
library(onemap)
library(ggplot2)

# load data 
F2.data <- read.mapmaker(file="F2_geno_for_one_map_final.txt")

# estimate two-point rf
twopts.f2.LOD3_rf0.5 <- rf.2pts(F2.data)

# save & export data 
save(twopts.f2.LOD3_rf0.5, file='twopts.LOD3_rf0.5.Rdata')





