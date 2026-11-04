#' Explore which features LymphGen selects for a gene panel or data type
#'
#' LymphGen doesn't have a fixed feature set: its `Model9()` refits the model
#' from the training data every time it runs, using only the genes and data
#' types available, and selects features for nine submodels (one per
#' combination of copy number and BCL2/BCL6 translocation status being
#' available or not, plus A53). Smaller gene lists or missing data types can
#' therefore bring in features that the complete model doesn't use (e.g. BCL2
#' Synon as a stand-in for the BCL2 translocation). This function reruns
#' LymphGen's own feature selection for one or more scenarios and reports
#' the features each submodel selects and how they differ from the
#' complete-data model.
#'
#' The training data are those of [load_lymphgen_training()] (not distributed with the package);
#' LymphGen's code (`Model9.R`) is sourced from your own copy of the LymphGen
#' R distribution, which isn't distributed with this package.
#'
#' @param lymphgen_dir Directory of the LymphGen R distribution, containing
#'   `Model9.R`.
#' @param scenarios Named list of scenarios, each a list with any of:
#'   * `mut_genes`: genes with mutation data (symbols, as in LymphGen's
#'     training annotation, or Entrez IDs). NULL (default): every gene.
#'   * `cn_genes`: genes with copy number data, as for `mut_genes`. NULL
#'     (default): every gene.
#'   * `include_arms`: whether arm-level copy number features are available
#'     (default TRUE; they're only used by the A53 model).
#'   * `cn_types`: which copy number calls are available: "all" (default),
#'     "homdel_amp" (only HOMDEL and AMP) or "hetloss_gain" (only HETLOSS
#'     and GAIN), LymphGen's CGH.class 0, 2 and 3.
#'   * `has_trunc`: whether truncating mutations are distinguished (default
#'     TRUE). Without them, LymphGen drops its TRUNC features.
#'   The default, `list(complete = list())`, is the complete-data model.
#' @param submodels Which submodels to fit: any of "Full", "NoCGH",
#'   "NoBCL2Fus", "NoBCL2CGH", "NoBCL6Fus", "NoBCL6CGH", "NoFus", "NoFusCGH"
#'   and "A53" (default all). Each scenario takes roughly 10-20 seconds per
#'   submodel.
#' @param lymphgen_data Path to the training data built by
#'   `data-raw/lymphgen_training.R` (see [load_lymphgen_training()]). Default:
#'   the `GAMBLR.predict.lymphgen_training` option, or else a copy in the
#'   package's `extdata`, if any.
#' @param verbose Print progress messages (default TRUE).
#'
#' @return A list with
#'   * `features`: one row per scenario, submodel and selected feature:
#'     `scenario`, `submodel`, `subtype`, `classifier_feature` (named as in
#'     GAMBLR.utils' classifier feature registry, e.g.
#'     "TNFAIP3:TRUNC|HOMDEL"), `gene`, `entrez_id`, `types` (the feature's
#'     levels, core first), and `in_complete_model` (whether the same
#'     submodel selects it with complete data).
#'   * `changes`: the features each scenario's submodels add or lose
#'     compared with the complete-data model: `scenario`, `submodel`,
#'     `subtype`, `classifier_feature`, `change` ("added" or "lost").
#'
#' @export
#'
#' @examples
#' \dontrun{
#' lyseq = readr::read_tsv(system.file("extdata", "lyseq_genes.tsv",
#'                                     package = "GAMBLR.predict"))
#' sel = lymphgen_feature_selection(
#'   "~/Downloads/Lymphgen for R 2",
#'   scenarios = list(lyseq = list(mut_genes = lyseq[[1]], cn_genes = lyseq[[1]]),
#'                    no_trunc = list(has_trunc = FALSE)),
#'   submodels = c("Full", "NoCGH"))
#' subset(sel$changes, scenario == "lyseq")
#' }
lymphgen_feature_selection = function(lymphgen_dir,
                                      scenarios = list(complete = list()),
                                      submodels = c("Full", "NoCGH", "NoBCL2Fus", "NoBCL2CGH",
                                                    "NoBCL6Fus", "NoBCL6CGH", "NoFus",
                                                    "NoFusCGH", "A53"),
                                      lymphgen_data = getOption("GAMBLR.predict.lymphgen_training", system.file("extdata", "lymphgen_training.rds",
                                                                  package = "GAMBLR.predict")),
                                      verbose = TRUE){
  all_submodels = c("Full", "NoCGH", "NoBCL2Fus", "NoBCL2CGH", "NoBCL6Fus",
                    "NoBCL6CGH", "NoFus", "NoFusCGH", "A53")
  unknown = setdiff(submodels, all_submodels)
  if(length(unknown) > 0){
    stop("unknown submodel: ", paste(unknown, collapse = ", "))
  }
  model_file = file.path(lymphgen_dir, "Model9.R")
  if(!file.exists(model_file)){
    stop("Model9.R not found in ", lymphgen_dir)
  }
  if(!nzchar(lymphgen_data) || !file.exists(lymphgen_data)){
    stop("LymphGen training data not found", if(nzchar(lymphgen_data)) paste0(": ", lymphgen_data) else "", ". ",
         "They aren't distributed with GAMBLR.predict: build them from the LymphGen R distribution with ",
         "data-raw/lymphgen_training.R (in the GAMBLR.predict repository) and pass the file as lymphgen_data, ",
         "or set options(GAMBLR.predict.lymphgen_training = <path>)")
  }
  if(is.null(names(scenarios)) || any(!nzchar(names(scenarios)))){
    stop("scenarios must be a named list")
  }
  allowed = c("mut_genes", "cn_genes", "include_arms", "cn_types", "has_trunc")
  for(nm in names(scenarios)){
    bad = setdiff(names(scenarios[[nm]]), allowed)
    if(length(bad) > 0){
      stop("scenario ", nm, ": unknown setting(s) ", paste(bad, collapse = ", "),
           ". Allowed: ", paste(allowed, collapse = ", "))
    }
  }

  lg = readRDS(lymphgen_data)
  uber = lg$uber; study = lg$study
  label = study$Stage3.genclass
  label[label == "SDT"] = "ST2"
  n_feat = lg$fullmat$dim[1]
  fullmat = matrix(0L, n_feat, lg$fullmat$dim[2])
  fullmat[lg$fullmat$ones] = 1L
  fullmat[lg$fullmat$na] = NA_integer_
  fullmat = as.data.frame(fullmat)

  # LymphGen's code uses globals and its default gene lists refer to an
  # undefined object, so run it in its own environment with explicit lists
  env = new.env()
  sys.source(model_file, envir = env)
  env$Full.Uber.2019 = uber
  env$Fullmat.2019 = fullmat

  model_features = function(models, scenario){
    out = NULL
    for(sm in names(models)){
      model = models[[sm]]
      if(is.null(model)) next
      mt = model$Mod.type
      keep = !is.na(mt$out) & mt$out > 0
      if(!any(keep)) next
      idx = as.matrix(mt[keep, intersect(c("out", "ot1", "ot1.1"), colnames(mt)), drop = FALSE])
      types = apply(idx, 1, function(r) paste(uber$Feature.type[r[!is.na(r) & r > 0]], collapse = "|"))
      gene = uber$Gene.Symbol[mt$out[keep]]
      out = rbind(out, data.frame(scenario = scenario, submodel = sm,
                                  subtype = model$lst[mt$i[keep]],
                                  classifier_feature = paste(gene, types, sep = ":"),
                                  gene = gene, entrez_id = uber$ENTREZ.ID[mt$out[keep]],
                                  types = types, stringsAsFactors = FALSE))
    }
    out
  }
  reference = model_features(lg$models[submodels], "complete")
  reference_key = paste(reference$submodel, reference$subtype, reference$classifier_feature)

  # genes given as symbols (old or new, e.g. HIST1H1D or H1-3; LymphGen
  # uses the older ones) or Entrez IDs, as LymphGen's Entrez IDs
  to_entrez = function(genes, what){
    if(is.numeric(genes)) return(genes)
    ids = lapply(genes, function(g){
      id = uber$ENTREZ.ID[match(GAMBLR.utils::expand_gene_aliases(g), uber$Gene.Symbol)]
      id[!is.na(id)]
    })
    missing = genes[lengths(ids) == 0]
    if(length(missing) > 0){
      warning(length(missing), " ", what, " not in LymphGen's training annotation: ",
              paste(utils::head(missing, 10), collapse = ", "),
              if(length(missing) > 10) ", ..." else "")
    }
    unique(unlist(ids))
  }
  all_genes = unique(uber$ENTREZ.ID[uber$ENTREZ.ID > 0])
  arms = unique(uber$ENTREZ.ID[uber$ENTREZ.ID < 0])

  features = NULL
  for(nm in names(scenarios)){
    s = scenarios[[nm]]
    mut = if(is.null(s$mut_genes)) all_genes else to_entrez(s$mut_genes, "mutation genes")
    cn = if(is.null(s$cn_genes)) all_genes else to_entrez(s$cn_genes, "copy number genes")
    # MYD88 L265P is its own pseudo-gene in LymphGen
    if(4615 %in% mut) mut = unique(c(mut, 1000000004))
    include_arms = if(is.null(s$include_arms)) TRUE else s$include_arms
    if(include_arms) cn = unique(c(cn, arms))
    cn_types = if(is.null(s$cn_types)) "all" else s$cn_types
    cgh_class = c(all = 0, homdel_amp = 2, hetloss_gain = 3)[cn_types]
    if(is.na(cgh_class)) stop("scenario ", nm, ": cn_types must be all, homdel_amp or hetloss_gain")
    has_trunc = if(is.null(s$has_trunc)) TRUE else s$has_trunc
    forcedrop = rep(FALSE, nrow(uber))
    if(!include_arms) forcedrop[uber$ENTREZ.ID < 0] = TRUE
    if(verbose){
      message("scenario ", nm, ": ", length(mut), " mutation genes, ", length(cn),
              " copy number genes/arms; fitting ", length(submodels), " submodel(s)...")
    }
    invisible(utils::capture.output(suppressWarnings(
      models <- env$Model9(lab = label, Full.Uber.2019 = uber, Fullmat.2019 = fullmat,
                           ForceCGH = c(6689, 5966, 1029, 4853, 7874),
                           CGH.class = unname(cgh_class), Has.CGH = study$CGH == 1,
                           modelcheck = match(submodels, all_submodels),
                           mut.genelist = mut, CGH.genelist = cn,
                           Has.Trunc = has_trunc, forcedrop = forcedrop))))
    names(models) = all_submodels
    features = rbind(features, model_features(models[submodels], nm))
  }
  features$in_complete_model = paste(features$submodel, features$subtype,
                                     features$classifier_feature) %in% reference_key

  changes = NULL
  for(nm in names(scenarios)){
    for(sm in submodels){
      these = features[features$scenario == nm & features$submodel == sm, ]
      ref = reference[reference$submodel == sm, ]
      added = these[!these$in_complete_model, c("subtype", "classifier_feature")]
      lost = ref[!paste(ref$subtype, ref$classifier_feature) %in%
                   paste(these$subtype, these$classifier_feature), c("subtype", "classifier_feature")]
      if(nrow(added) + nrow(lost) == 0) next
      changes = rbind(changes, data.frame(scenario = nm, submodel = sm,
                                          rbind(added, lost),
                                          change = rep(c("added", "lost"), c(nrow(added), nrow(lost))),
                                          stringsAsFactors = FALSE))
    }
  }
  if(is.null(changes)){
    changes = data.frame(scenario = character(), submodel = character(), subtype = character(),
                         classifier_feature = character(), change = character())
  }
  rownames(features) = NULL; rownames(changes) = NULL
  list(features = features, changes = changes)
}
