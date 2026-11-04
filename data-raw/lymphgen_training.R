# Builds inst/extdata/lymphgen_training.rds, used by load_lymphgen_training(),
# from the LymphGen R distribution (Wright et al. 2020, Cancer Cell; "Lymphgen
# for R": Model9.R, Fullmat.2019.txt, Full.Uber.2019.txt, Study.sampR.txt).
#
# Contents:
#   uber: feature annotation (Full.Uber.2019.txt; Model9 needs all columns)
#   study: training samples (Study.sampR.txt)
#   fullmat: the training matrix (Fullmat.2019.txt, 163,427 features x 574
#     samples, 0/1/NA) as sorted linear cell indices of its 1s and NAs
#   models: LymphGen's models from its own Model9() on the training data
#     with every gene available, all nine submodels, named Full, NoCGH,
#     NoBCL2Fus, NoBCL2CGH, NoBCL6Fus, NoBCL6CGH, NoFus, NoFusCGH and A53
#   prefilter: TRUE for features passing Model9's prefilters (Get.Featureset)
#     in either model, with every gene available
#
# Usage: Rscript data-raw/lymphgen_training.R "~/Downloads/Lymphgen for R 2"

args = commandArgs(trailingOnly = TRUE)
lymphgen_dir = if(length(args)) args[1] else "~/Downloads/Lymphgen for R 2"
rd = function(f) read.table(file.path(lymphgen_dir, f), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
uber = rd("Full.Uber.2019.txt")
fullmat = rd("Fullmat.2019.txt")
study = rd("Study.sampR.txt")
label = study$Stage3.genclass
label[label == "SDT"] = "ST2"

# Model9.R uses globals and its default gene lists refer to an undefined
# Full.Uber, so pass the gene lists and run it in its own environment
env = new.env()
sys.source(file.path(lymphgen_dir, "Model9.R"), envir = env)
env$Full.Uber.2019 = uber; env$Fullmat.2019 = fullmat
genes = unique(uber[, 3])
invisible(utils::capture.output(
  models <- env$Model9(lab = label, Full.Uber.2019 = uber, Fullmat.2019 = fullmat,
                       Has.CGH = study$CGH == 1, modelcheck = 1:9,
                       mut.genelist = genes, CGH.genelist = genes)))
names(models) = c("Full", "NoCGH", "NoBCL2Fus", "NoBCL2CGH", "NoBCL6Fus",
                  "NoBCL6CGH", "NoFus", "NoFusCGH", "A53")

# the prefilters Model9 applies with every gene available
force = c(6689, 5966, 1029, 4853, 7874)
no_drop = rep(FALSE, nrow(uber))
keep_full = env$Get.Featureset(uber, ForceCGH = force, drop.CGH = 0, custom.Drop = no_drop)
keep_a53 = env$Get.Featureset(uber, exclude.Other = FALSE, ForceCGH = force, Stage = 2, custom.Drop = no_drop)
prefilter = rowSums(keep_full) > 0 | rowSums(keep_a53) > 0

m = as.matrix(fullmat)
lymphgen_training = list(
  uber = uber,
  study = study,
  fullmat = list(dim = dim(m), ones = which(m == 1), na = which(is.na(m))),
  models = models,
  prefilter = prefilter
)
saveRDS(lymphgen_training, "inst/extdata/lymphgen_training.rds", compress = "xz")
