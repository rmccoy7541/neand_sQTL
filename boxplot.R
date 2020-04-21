library(data.table)

# box plot script

wholeblood <- fread("Whole_Blood_NL_isos.txt")





library("mlmm.gwas")
data("mlmm.gwas.AD")

XX = list(Xa)
KK = list(K.add)

# GWAS
res_mlmm <- mlmm_allmodels(floweringDateAD, XX, KK)
manhattan.plot(res_mlmm)

# Model selection
sel_XX <- frommlmm_toebic(XX, res_mlmm)
res.eBIC <- eBIC_allmodels(floweringDateAD, sel_XX, KK, ncol(Xa))

# Effects estimations with the selected model
sel_XXclass <- fromeBICtoEstimation(XX, res.eBIC)
eff.estimations <- Estimation_allmodels(floweringDateAD, sel_XXclass, KK)
genotypes.boxplot(Xa, floweringDateAD, effects = eff.estimations)
# }