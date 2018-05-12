library(tcltk)
library(tkrplot)
library(onemap)  
library(ggplot2)
library(dplyr)

F2.data <- read.mapmaker(file="F2_geno_for_one_map_final.txt")
load("twopts.LOD3_rf0.5.Rdata")
load("group.AC.Rdata")

# order marker within each chromosome

mark.f2 <- list()
LG.f2.ord <- list()

for (i in 1:length(group.AC)) {
	mark.f2[[i]] <- make.seq(twopts.f2.LOD3_rf0.5, group.AC[[i]])
	LG.f2.ord[[i]] <- order.seq(input.seq = mark.f2[[i]], n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 3,
                        draw.try = T, wait = 1,
                        touchdown = T)
}

save(mark.f2, LG.f2.ord, file = "LG.f2.ord.AC.Rdata")




