library(qtl)

LG.f2.after.crossover <- read.cross("csvsr", genfile = "/share/malooflab/Ruijuan/F2/data/LG.f2.madmapper.final_gen_revised_flipped.csv",
                     phefile = "/share/malooflab/Ruijuan/F2/data/growth_model_trait.2.csv",
                     genotypes = c("AA", "AB", "BB")) # the only problem is that 44 phenotypes were read instead of 43, need to figure out why later

LG.f2.after.crossover <- sim.geno(LG.f2.after.crossover,step=1,n.draws=32)
LG.f2.after.crossover <- calc.genoprob(LG.f2.after.crossover,step=1)

# run scantwo for all traits at once 
system.time(
scantwo.imp <- 
lapply(seq_along(LG.f2.after.crossover$pheno)[1:12], function(trait) {
  print(trait)
  scantwo(LG.f2.after.crossover,pheno.col=trait,method="hk")
}) 
)

names(scantwo.imp) <- colnames(LG.f2.after.crossover$pheno)[1:12]

# save output 
save(scantwo.imp, file = "/share/malooflab/Ruijuan/F2/QTL_analysis/output/scantwo.imp.v2.growth.model.flipped.Rdata")
