#' @title DLBCLone class evidence weights
#'
#' @description Derive, for each class, the features that support it and a
#' weight for each, from training features and labels. Used by the DLBCLone
#' evidence gate ([DLBCLone_evidence_gate]), which keeps a class call only when
#' the sample carries enough weighted evidence for that class.
#'
#' @details A feature supports a class when it is associated with the class in
#' the training samples: a one-vs-rest Fisher test on presence (value > 0)
#' against the class's samples, Benjamini-Hochberg adjusted over all features
#' and classes, with `q < max_q` and an odds ratio above 1. Its weight is
#' `strong_weight` (default 2) when the feature is specific to the class, i.e.
#' at least `strong_ppv` of the training samples carrying it are in the class
#' (a defining lesion such as MYD88 L265P for MCD), and `weight` (default 1)
#' otherwise. A feature can support more than one class. Specificity is used
#' rather than the odds ratio because features the labelling classifier itself
#' uses (e.g. BCL6 translocations for LymphGen's BN2) can reach a high odds
#' ratio while being common outside the class.
#'
#' @param features Data frame or matrix of features (rows = samples, sample IDs
#' as row names), as given to [DLBCLone_optimize_params]. Meta-feature columns
#' ending in "_feats" (from `core_features`) are ignored.
#' @param labels Named character vector of class labels (names = sample IDs),
#' or a data frame with `sample_id` and `truth_column`.
#' @param classes Classes to derive evidence for (default: every label except
#' `other_class`).
#' @param truth_column With a data frame `labels`: the column holding the labels.
#' @param other_class The label meaning "no class" (default "Other"); never a
#' class here.
#' @param max_q Benjamini-Hochberg q-value below which a feature supports a
#' class (default 0.05).
#' @param strong_ppv Share of the feature's carriers that must be in the class
#' for the feature to get `strong_weight` (default 0.7).
#' @param weight,strong_weight The two weights (defaults 1 and 2).
#'
#' @return A data frame with one row per supporting feature and class:
#' `class`, `feature`, `weight`, `n` (training samples carrying it), `ppv`,
#' `odds_ratio`, `q`.
#'
#' @examples
#' \dontrun{
#' ev = DLBCLone_class_evidence(features, metadata, truth_column = "lymphgen")
#' }
#' @import dplyr
#' @export
DLBCLone_class_evidence = function(features,
                                   labels,
                                   classes = NULL,
                                   truth_column = "lymphgen",
                                   other_class = "Other",
                                   max_q = 0.05,
                                   strong_ppv = 0.7,
                                   weight = 1,
                                   strong_weight = 2){
  features = .dlbclone_feature_rows(features)
  if(is.data.frame(labels)){
    if(!all(c("sample_id", truth_column) %in% colnames(labels))){
      stop("labels needs sample_id and ", truth_column, " columns")
    }
    labels = stats::setNames(as.character(labels[[truth_column]]), labels$sample_id)
  }
  ids = intersect(rownames(features), names(labels))
  if(length(ids) == 0) stop("no samples shared by features (row names) and labels")
  y = labels[ids]
  keep = !is.na(y)
  ids = ids[keep]; y = y[keep]
  x = features[ids, !grepl("_feats$", colnames(features)), drop = FALSE] > 0
  if(is.null(classes)) classes = setdiff(sort(unique(y)), other_class)
  classes = setdiff(classes, other_class)
  tests = lapply(classes, function(k){
    in_k = y == k
    do.call(rbind, lapply(colnames(x), function(f){
      hit = x[, f]
      if(!any(hit) || all(hit)) return(NULL)
      ft = stats::fisher.test(table(factor(in_k, c(TRUE, FALSE)), factor(hit, c(TRUE, FALSE))))
      data.frame(class = k, feature = f, n = sum(hit), ppv = mean(in_k[hit]),
                 odds_ratio = unname(ft$estimate), p = ft$p.value, stringsAsFactors = FALSE)
    }))
  })
  tests = do.call(rbind, tests)
  if(is.null(tests) || nrow(tests) == 0){
    return(data.frame(class = character(), feature = character(), weight = numeric(), n = integer(),
                      ppv = numeric(), odds_ratio = numeric(), q = numeric()))
  }
  tests$q = stats::p.adjust(tests$p, "BH")
  out = tests[tests$q < max_q & tests$odds_ratio > 1, , drop = FALSE]
  out$weight = ifelse(out$ppv >= strong_ppv, strong_weight, weight)
  out = out[order(match(out$class, classes), -out$weight, out$q), c("class", "feature", "weight", "n", "ppv", "odds_ratio", "q")]
  rownames(out) = NULL
  out
}

#' @title DLBCLone evidence gate
#'
#' @description Keep a DLBCLone class call only when the sample carries enough
#' weighted evidence for that class; otherwise call it `other_class`.
#'
#' @details A sample's evidence for a class is the sum of the weights
#' ([DLBCLone_class_evidence]) of the class's supporting features the sample
#' carries (value > 0). A call of class X is kept when the evidence for X is at
#' least `min_evidence` (a single number, or a named vector giving each class
#' its own minimum; classes not named use the "default" element if given,
#' else 2). With the default weights and minimum of 2, a call needs two
#' supporting features or one class-specific (weight 2) feature. Calls of
#' classes without evidence weights (e.g. composites) are kept unchanged and
#' get NA evidence. The input call column is left as it is; the gate adds
#' columns.
#'
#' @param predictions Data frame of predictions with `sample_id` and
#' `pred_column`.
#' @param features Data frame or matrix of the samples' features (row names =
#' sample IDs), unweighted or weighted (only presence counts).
#' @param class_evidence Output of [DLBCLone_class_evidence] (or a data frame
#' with `class`, `feature` and `weight`).
#' @param min_evidence Minimum weighted evidence to keep a call (default 2);
#' see details.
#' @param pred_column Column holding the calls to gate (default
#' "DLBCLone_wo").
#' @param out_column Column for the gated calls (default
#' `paste0(pred_column, "g")`, e.g. "DLBCLone_wog").
#' @param other_class The label meaning "no class" (default "Other").
#'
#' @return `predictions` with three columns added: `<out_column>` (the gated
#' call), `DLBCLone_evidence` (the evidence for the class in `pred_column`)
#' and `DLBCLone_evidence_min` (the minimum it was held to).
#'
#' @examples
#' \dontrun{
#' ev = DLBCLone_class_evidence(features, metadata, truth_column = "lymphgen")
#' gated = DLBCLone_evidence_gate(model$predictions, features, ev, min_evidence = c(default = 2, ST2 = 3))
#' }
#' @export
DLBCLone_evidence_gate = function(predictions,
                                  features,
                                  class_evidence,
                                  min_evidence = 2,
                                  pred_column = "DLBCLone_wo",
                                  out_column = NULL,
                                  other_class = "Other"){
  if(is.null(out_column)) out_column = paste0(pred_column, "g")
  if(!all(c("sample_id", pred_column) %in% colnames(predictions))){
    stop("predictions needs sample_id and ", pred_column, " columns")
  }
  if(!all(c("class", "feature", "weight") %in% colnames(class_evidence))){
    stop("class_evidence needs class, feature and weight columns")
  }
  features = .dlbclone_feature_rows(features)
  calls = as.character(predictions[[pred_column]])
  ids = as.character(predictions$sample_id)
  missing_ids = setdiff(ids[calls %in% class_evidence$class], rownames(features))
  if(length(missing_ids)){
    stop(length(missing_ids), " sample(s) with a class call have no row in features, e.g. ", missing_ids[1])
  }
  missing_features = setdiff(class_evidence$feature, colnames(features))
  if(length(missing_features)){
    message("evidence gate: ", length(missing_features), " supporting feature(s) not in features, counted as absent: ",
            paste(utils::head(missing_features, 10), collapse = ", "))
  }
  min_for = function(k){
    if(length(min_evidence) == 1 && is.null(names(min_evidence))) return(unname(min_evidence))
    if(k %in% names(min_evidence)) return(unname(min_evidence[[k]]))
    if("default" %in% names(min_evidence)) return(unname(min_evidence[["default"]]))
    2
  }
  evidence = rep(NA_real_, length(ids))
  threshold = rep(NA_real_, length(ids))
  for(k in unique(calls[calls %in% class_evidence$class])){
    ev = class_evidence[class_evidence$class == k & class_evidence$feature %in% colnames(features), , drop = FALSE]
    rows = which(calls == k)
    present = as.matrix(features[ids[rows], ev$feature, drop = FALSE]) > 0
    present[is.na(present)] = FALSE
    evidence[rows] = as.vector(present %*% ev$weight)
    threshold[rows] = min_for(k)
  }
  gated = ifelse(!is.na(evidence) & evidence < threshold, other_class, calls)
  predictions[[out_column]] = gated
  predictions$DLBCLone_evidence = evidence
  predictions$DLBCLone_evidence_min = threshold
  n_gated = sum(gated != calls, na.rm = TRUE)
  message("evidence gate: ", n_gated, " of ", sum(calls %in% class_evidence$class), " class calls set to ",
          other_class, " (", out_column, ")")
  predictions
}

# features with sample IDs as row names (a sample_id column, if present, becomes the row names)
.dlbclone_feature_rows = function(features){
  features = as.data.frame(features)
  if("sample_id" %in% colnames(features)){
    rownames(features) = features$sample_id
    features$sample_id = NULL
  }
  features
}
