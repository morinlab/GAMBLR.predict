#' @title DLBCLone cross-validation folds
#'
#' @description Assign samples to cross-validation folds, stratified by class
#' and keeping each group (e.g. a patient's tumours) in one fold.
#'
#' @details Groups are dealt to folds class by class in a random order, so each
#' class is spread as evenly as possible across folds. A group's class is the
#' label of its first sample.
#'
#' @param metadata Data frame with `sample_id` and `truth_column`.
#' @param truth_column Column holding the class labels (default "lymphgen").
#' @param k Number of folds (default 5).
#' @param group_column Optional column whose values must stay together (e.g.
#' "patient_id"); default NULL: each sample is its own group.
#' @param seed Random seed (default 12345).
#'
#' @return Data frame with `sample_id` and `fold` (1 to `k`).
#'
#' @examples
#' md = data.frame(sample_id = paste0("s", 1:20), lymphgen = rep(c("EZB", "MCD"), 10))
#' DLBCLone_cv_folds(md, k = 5)
#' @export
DLBCLone_cv_folds = function(metadata,
                             truth_column = "lymphgen",
                             k = 5,
                             group_column = NULL,
                             seed = 12345){
  md = as.data.frame(metadata)
  if(!all(c("sample_id", truth_column) %in% colnames(md))) stop("metadata needs sample_id and ", truth_column)
  if(anyDuplicated(md$sample_id)) stop("metadata has duplicated sample_id values")
  if(!is.null(group_column) && !group_column %in% colnames(md)) stop("no ", group_column, " column in metadata")
  grp = if(is.null(group_column)) md$sample_id else as.character(md[[group_column]])
  groups = unique(grp)
  label = as.character(md[[truth_column]])[match(groups, grp)]
  label[is.na(label)] = "NA"
  set.seed(seed)
  fold = integer(length(groups))
  dealt = 0
  for(l in sort(unique(label))){
    idx = which(label == l)
    idx = idx[sample.int(length(idx))]
    fold[idx] = ((dealt + seq_along(idx) - 1) %% k) + 1
    dealt = dealt + length(idx)
  }
  data.frame(sample_id = md$sample_id, fold = fold[match(grp, groups)], stringsAsFactors = FALSE)
}

#' @title DLBCLone cross-validation
#'
#' @description Estimate how DLBCLone performs on samples it wasn't trained on:
#' for each fold, train a model on the other folds ([make_and_annotate_umap],
#' [DLBCLone_optimize_params]) and classify the held-out fold with
#' [DLBCLone_predict].
#'
#' @details Everything learned from labels happens inside each fold, on its
#' training samples only: the UMAP, the k/threshold optimisation, the evidence
#' gate's weights (see [DLBCLone_class_evidence]) and, with
#' `select_features = TRUE`, the feature selection. Feature selection keeps the
#' features associated with a class (those [DLBCLone_class_evidence] returns,
#' with `select_args`); features absent from a fold's training samples are
#' dropped. Samples whose label isn't in `truth_classes` (e.g. composite
#' calls) are never trained on but are classified when their fold is held out.
#'
#' Folds can be run separately (`run_folds`), e.g. one fold per R process to
#' keep memory down, and the results combined with [DLBCLone_cv_combine].
#'
#' @param features Data frame of features (rows = samples, sample IDs as row
#' names or a `sample_id` column).
#' @param metadata Data frame with `sample_id` and `truth_column` (and
#' `group_column`, if used).
#' @param truth_classes Classes to train on (including `other_class` if
#' Other is a training label).
#' @param truth_column Column holding the labels (default "lymphgen").
#' @param other_class The label meaning "no class" (default "Other").
#' @param k Number of folds (default 5); ignored when `folds` is given.
#' @param folds Optional fold assignment (data frame with `sample_id` and
#' `fold`, as [DLBCLone_cv_folds] returns); default: [DLBCLone_cv_folds] with
#' `k`, `group_column` and `seed`. Use the same folds to compare models.
#' @param group_column Optional column of `metadata` whose values stay in one
#' fold (e.g. "patient_id").
#' @param run_folds Which folds to run (default: all).
#' @param select_features If TRUE, select features within each fold (see
#' details). Default FALSE.
#' @param select_args List of arguments to [DLBCLone_class_evidence] for the
#' feature selection (e.g. `list(max_q = 0.01)`).
#' @param umap_args,optimize_args,predict_args Lists of further arguments to
#' [make_and_annotate_umap], [DLBCLone_optimize_params] (e.g.
#' `list(optimize_for_other = TRUE)`, `evidence_min`) and [DLBCLone_predict].
#' @param seed Random seed for the folds and training (default 12345).
#'
#' @return A list: `predictions` (one row per held-out sample: `sample_id`,
#' `fold`, the label as `truth`, and DLBCLone's calls and scores, including
#' `DLBCLone_wo` and, with the evidence gate, `DLBCLone_wog`; held-out samples
#' without any of the fold's features are called `other_class`, flagged in
#' `no_features`), `params` (each
#' fold's chosen parameters and number of features), `features_used` (per
#' fold), `folds`, and `settings`.
#'
#' @examples
#' \dontrun{
#' cv = DLBCLone_cross_validate(features, metadata, truth_classes = c("EZB", "MCD", "ST2", "BN2", "Other"),
#'                              group_column = "patient_id", select_features = TRUE,
#'                              optimize_args = list(optimize_for_other = TRUE))
#' DLBCLone_cv_summary(cv)
#' }
#' @export
DLBCLone_cross_validate = function(features,
                                   metadata,
                                   truth_classes,
                                   truth_column = "lymphgen",
                                   other_class = "Other",
                                   k = 5,
                                   folds = NULL,
                                   group_column = NULL,
                                   run_folds = NULL,
                                   select_features = FALSE,
                                   select_args = list(),
                                   umap_args = list(),
                                   optimize_args = list(),
                                   predict_args = list(),
                                   seed = 12345){
  features = .dlbclone_feature_rows(features)
  md = as.data.frame(metadata)
  md = md[md$sample_id %in% rownames(features), , drop = FALSE]
  if(nrow(md) == 0) stop("no samples shared by features and metadata")
  if(is.null(folds)) folds = DLBCLone_cv_folds(md, truth_column = truth_column, k = k, group_column = group_column, seed = seed)
  folds = folds[folds$sample_id %in% md$sample_id, , drop = FALSE]
  if(is.null(run_folds)) run_folds = sort(unique(folds$fold))
  trainable = md$sample_id[md[[truth_column]] %in% truth_classes]
  preds = list(); params = list(); used = list()
  for(f in run_folds){
    test_ids = folds$sample_id[folds$fold == f]
    train_md = md[md$sample_id %in% setdiff(trainable, test_ids), , drop = FALSE]
    message("== fold ", f, ": training on ", nrow(train_md), " samples, holding out ", length(test_ids))
    x_train = features[train_md$sample_id, , drop = FALSE]
    if(select_features){
      ev = do.call(DLBCLone_class_evidence,
                   c(list(features = x_train, labels = train_md, classes = setdiff(truth_classes, other_class),
                          truth_column = truth_column, other_class = other_class), select_args))
      keep = union(unique(ev$feature), grep("_feats$", colnames(x_train), value = TRUE))
      x_train = x_train[, colnames(x_train) %in% keep, drop = FALSE]
    }
    x_train = x_train[, colSums(x_train > 0, na.rm = TRUE) > 0, drop = FALSE]
    message("fold ", f, ": ", ncol(x_train), " features")
    umap_out = do.call(make_and_annotate_umap,
                       c(list(df = x_train, metadata = train_md, truth_column = truth_column, seed = seed), umap_args))
    opt = do.call(DLBCLone_optimize_params,
                  c(list(combined_mutation_status_df = x_train, metadata_df = train_md, umap_out = umap_out,
                         truth_classes = truth_classes, truth_column = truth_column, other_class = other_class,
                         seed = seed), optimize_args))
    x_test = features[test_ids, colnames(x_train), drop = FALSE]
    pr = do.call(DLBCLone_predict, c(list(mutation_status = x_test, optimized_model = opt, seed = seed), predict_args))$prediction
    pr = as.data.frame(pr)
    # held-out samples DLBCLone_predict leaves out (no features) are Other
    dropped = setdiff(test_ids, pr$sample_id)
    pr$no_features = FALSE
    if(length(dropped)){
      filler = data.frame(sample_id = dropped, no_features = TRUE, stringsAsFactors = FALSE)
      for(v in intersect(c("DLBCLone_w", "DLBCLone_wo", "DLBCLone_wog"), colnames(pr))) filler[[v]] = other_class
      pr = dplyr::bind_rows(pr, filler)
      message("fold ", f, ": ", length(dropped), " held-out sample(s) without features called ", other_class)
    }
    pr$fold = f
    pr$truth = as.character(md[[truth_column]])[match(pr$sample_id, md$sample_id)]
    lead = c("sample_id", "fold", "truth")
    preds[[as.character(f)]] = pr[, c(lead, setdiff(colnames(pr), lead)), drop = FALSE]
    bp = as.data.frame(opt$best_params)
    bp$fold = f; bp$n_features = ncol(x_train); bp$n_train = nrow(train_md); bp$n_test = length(test_ids)
    params[[as.character(f)]] = bp
    used[[as.character(f)]] = colnames(x_train)
    rm(opt, umap_out); invisible(gc())
  }
  list(predictions = dplyr::bind_rows(preds), params = dplyr::bind_rows(params), features_used = used,
       folds = folds, settings = list(truth_classes = truth_classes, truth_column = truth_column,
                                      other_class = other_class, select_features = select_features, seed = seed))
}

#' @title Combine DLBCLone cross-validation runs
#'
#' @description Combine [DLBCLone_cross_validate] results from separate runs
#' (e.g. one fold each) into one.
#'
#' @param cv_list List of [DLBCLone_cross_validate] results.
#'
#' @return One result in the same form, with every fold's predictions and
#' parameters. Stops if a fold appears twice or the runs used different folds.
#' @export
DLBCLone_cv_combine = function(cv_list){
  folds = cv_list[[1]]$folds
  for(cv in cv_list[-1]){
    if(!isTRUE(all.equal(cv$folds[order(cv$folds$sample_id), ], folds[order(folds$sample_id), ], check.attributes = FALSE))){
      stop("the runs used different fold assignments")
    }
  }
  run = unlist(lapply(cv_list, function(cv) unique(cv$params$fold)))
  if(anyDuplicated(run)) stop("fold(s) present more than once: ", paste(unique(run[duplicated(run)]), collapse = ", "))
  list(predictions = dplyr::bind_rows(lapply(cv_list, `[[`, "predictions")),
       params = dplyr::bind_rows(lapply(cv_list, `[[`, "params")),
       features_used = do.call(c, lapply(cv_list, `[[`, "features_used")),
       folds = folds, settings = cv_list[[1]]$settings)
}

#' @title Summarise DLBCLone cross-validation
#'
#' @description Held-out performance per call column: agreement with the
#' labels, share classified, and per-class recall, precision and F1.
#'
#' @details Scored on the held-out samples whose label is in the training
#' classes; samples with other labels (e.g. composites) are left out. The
#' other class is not scored as a class: a sample labelled `other_class` but
#' given a class lowers that class's precision only.
#'
#' @param cv A [DLBCLone_cross_validate] (or [DLBCLone_cv_combine]) result.
#' @param pred_columns Call columns to score (default "DLBCLone_wo" and, if
#' present, "DLBCLone_wog").
#'
#' @return A list: `overall` (per column: samples, agreement, share
#' classified, labelled-Other samples given a class) and `per_class` (per
#' column and class: labelled samples, calls, recall, precision, F1).
#' @export
DLBCLone_cv_summary = function(cv, pred_columns = NULL){
  p = cv$predictions
  classes = cv$settings$truth_classes
  other = cv$settings$other_class
  if(is.null(pred_columns)) pred_columns = intersect(c("DLBCLone_wo", "DLBCLone_wog"), colnames(p))
  p = p[p$truth %in% classes, , drop = FALSE]
  real = setdiff(classes, other)
  overall = do.call(rbind, lapply(pred_columns, function(v){
    data.frame(call = v, samples = nrow(p), agreement = mean(p[[v]] == p$truth), classified = mean(p[[v]] != other),
               other_given_class = sum(p$truth == other & p[[v]] != other), stringsAsFactors = FALSE)
  }))
  per_class = do.call(rbind, lapply(pred_columns, function(v){
    do.call(rbind, lapply(real, function(kk){
      tp = sum(p$truth == kk & p[[v]] == kk); fp = sum(p$truth != kk & p[[v]] == kk); fn = sum(p$truth == kk & p[[v]] != kk)
      data.frame(call = v, class = kk, labelled = tp + fn, calls = tp + fp,
                 recall = if(tp + fn > 0) tp / (tp + fn) else NA_real_,
                 precision = if(tp + fp > 0) tp / (tp + fp) else NA_real_,
                 F1 = if(tp > 0) 2 * tp / (2 * tp + fp + fn) else 0, stringsAsFactors = FALSE)
    }))
  }))
  list(overall = overall, per_class = per_class)
}
