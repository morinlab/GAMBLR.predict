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

#' @title Recommend a k for DLBCLone
#'
#' @description Choose DLBCLone's k from held-out performance: cross-validate
#' each k in `ks` on the same folds and recommend one from the flat top of the
#' resulting curve. Train the final model with `min_k = max_k = ` the
#' recommended k.
#'
#' @details Within each fold (as [DLBCLone_cross_validate]: optional feature
#' selection on the fold's training samples, then UMAP), the UMAP is fitted
#' and the samples projected once, and every k reuses them: for each k,
#' [DLBCLone_optimize_params] (with `min_k = max_k = k`) tunes the
#' remaining settings on the training samples and [DLBCLone_predict]
#' classifies the held-out ones. The held-out calls of all folds are pooled
#' per k and scored on the samples whose label is a training class.
#'
#' The recommendation is the middle of the k values whose score is within
#' `tolerance` of the best (the median, rounded to a k in `ks`, the larger on
#' a tie): the curve is usually flat near its top, so its single best point
#' moves with small changes in the data. With `metric` "balanced_accuracy",
#' the score is the mean per-class balanced accuracy over all training
#' classes including `other_class` (what [DLBCLone_optimize_params]
#' maximises); with "agreement", the share of samples given their label.
#' Scores use `pred_column` (default the evidence-gated call when present).
#'
#' A k larger than the smallest class has few training samples per vote; a
#' warning is given when the recommended k exceeds the smallest class's
#' training samples in any fold.
#'
#' @inheritParams DLBCLone_cross_validate
#' @param ks The k values to compare (default 5 to 25, step 2).
#' @param metric "balanced_accuracy" (default) or "agreement".
#' @param tolerance How far below the best score a k may be and still count
#' as on the top of the curve (default 0.005).
#' @param pred_column Call column to score (default "DLBCLone_wog" when the
#' evidence gate is on, else "DLBCLone_wo").
#'
#' @return A list: `recommended_k`, `scores` (one row per k: the metric,
#' agreement, share classified, labelled-Other samples given a class, and
#' per-class F1), `near_best` (the k values within `tolerance`),
#' `predictions` (held-out calls, all k), `folds`, `plot` (score by k) and
#' `settings`.
#'
#' @examples
#' \dontrun{
#' rk = DLBCLone_recommend_k(features, metadata, truth_classes = c("EZB", "MCD", "ST2", "BN2", "Other"),
#'                           group_column = "patient_id", select_features = TRUE,
#'                           optimize_args = list(optimize_for_other = TRUE))
#' rk$recommended_k; rk$plot
#' }
#' @import ggplot2
#' @export
DLBCLone_recommend_k = function(features,
                                metadata,
                                truth_classes,
                                truth_column = "lymphgen",
                                other_class = "Other",
                                ks = seq(5, 25, 2),
                                k = 5,
                                folds = NULL,
                                group_column = NULL,
                                select_features = FALSE,
                                select_args = list(),
                                umap_args = list(),
                                optimize_args = list(),
                                predict_args = list(),
                                metric = "balanced_accuracy",
                                tolerance = 0.005,
                                pred_column = NULL,
                                seed = 12345){
  metric = match.arg(metric, c("balanced_accuracy", "agreement"))
  ks = sort(unique(ks))
  features = .dlbclone_feature_rows(features)
  md = as.data.frame(metadata)
  md = md[md$sample_id %in% rownames(features), , drop = FALSE]
  if(nrow(md) == 0) stop("no samples shared by features and metadata")
  if(is.null(folds)) folds = DLBCLone_cv_folds(md, truth_column = truth_column, k = k, group_column = group_column, seed = seed)
  folds = folds[folds$sample_id %in% md$sample_id, , drop = FALSE]
  trainable = md$sample_id[md[[truth_column]] %in% truth_classes]
  real = setdiff(truth_classes, other_class)
  preds = list()
  smallest = Inf
  for(f in sort(unique(folds$fold))){
    test_ids = folds$sample_id[folds$fold == f]
    train_md = md[md$sample_id %in% setdiff(trainable, test_ids), , drop = FALSE]
    smallest = min(smallest, min(table(factor(train_md[[truth_column]], real))))
    x_train = features[train_md$sample_id, , drop = FALSE]
    if(select_features){
      ev = do.call(DLBCLone_class_evidence,
                   c(list(features = x_train, labels = train_md, classes = real,
                          truth_column = truth_column, other_class = other_class), select_args))
      x_train = x_train[, colnames(x_train) %in% union(unique(ev$feature), grep("_feats$", colnames(x_train), value = TRUE)), drop = FALSE]
    }
    x_train = x_train[, colSums(x_train > 0, na.rm = TRUE) > 0, drop = FALSE]
    message("== fold ", f, ": ", nrow(train_md), " training samples, ", ncol(x_train), " features; k ", paste(ks, collapse = ","))
    umap_out = do.call(make_and_annotate_umap,
                       c(list(df = x_train, metadata = train_md, truth_column = truth_column, seed = seed), umap_args))
    # the training projection DLBCLone_optimize_params would make from umap_out, and the held-out projection
    # DLBCLone_predict would make: computed once per fold, reused for every k
    proj_train = make_and_annotate_umap(df = umap_out$features, umap_out = umap_out, min_dist = 0, n_neighbors = 55,
                                        n_epochs = 1500, seed = seed, metadata = train_md, ret_model = TRUE,
                                        metric = "cosine", join_column = "sample_id", na_vals = "drop")
    if("sample_metadata_no_features" %in% names(umap_out)) proj_train$sample_metadata_no_features = umap_out$sample_metadata_no_features
    x_test = features[test_ids, colnames(x_train), drop = FALSE]
    proj_test = make_and_annotate_umap(df = tibble::rownames_to_column(x_test, "sample_id"), umap_out = umap_out,
                                       ret_model = FALSE, seed = seed, join_column = "sample_id", na_vals = "drop")
    for(kk in ks){
      opt = do.call(DLBCLone_optimize_params,
                    c(list(combined_mutation_status_df = x_train, metadata_df = train_md, umap_out = umap_out,
                           umap_projection = proj_train, truth_classes = truth_classes, truth_column = truth_column,
                           other_class = other_class, seed = seed, min_k = kk, max_k = kk),
                      optimize_args[setdiff(names(optimize_args), c("min_k", "max_k"))]))
      opt$projection_train = proj_train$df
      pr = as.data.frame(do.call(DLBCLone_predict, c(list(mutation_status = x_test, optimized_model = opt, seed = seed,
                                                          test_projection = proj_test), predict_args))$prediction)
      dropped = setdiff(test_ids, pr$sample_id)
      pr$no_features = FALSE
      if(length(dropped)){
        filler = data.frame(sample_id = dropped, no_features = TRUE, stringsAsFactors = FALSE)
        for(v in intersect(c("DLBCLone_w", "DLBCLone_wo", "DLBCLone_wog"), colnames(pr))) filler[[v]] = other_class
        pr = dplyr::bind_rows(pr, filler)
      }
      keep = intersect(c("sample_id", "DLBCLone_w", "DLBCLone_wo", "DLBCLone_wog", "DLBCLone_evidence", "no_features"), colnames(pr))
      pr = pr[, keep, drop = FALSE]
      pr$k = kk; pr$fold = f
      pr$truth = as.character(md[[truth_column]])[match(pr$sample_id, md$sample_id)]
      pr$purity_DLBCLone_w = if(is.null(opt$best_params$purity_DLBCLone_w)) NA_real_ else opt$best_params$purity_DLBCLone_w
      preds[[paste(f, kk)]] = pr
      rm(opt); invisible(gc())
    }
    rm(umap_out, proj_train, proj_test); invisible(gc())
  }
  preds = dplyr::bind_rows(preds)
  if(is.null(pred_column)) pred_column = if("DLBCLone_wog" %in% colnames(preds)) "DLBCLone_wog" else "DLBCLone_wo"
  scored = preds[preds$truth %in% truth_classes, , drop = FALSE]
  scores = do.call(rbind, lapply(ks, function(kk){
    p = scored[scored$k == kk, , drop = FALSE]
    call = p[[pred_column]]
    row = data.frame(k = kk,
                     balanced_accuracy = .purity_accuracy(p$truth, call, truth_classes, "balanced_accuracy"),
                     agreement = mean(call == p$truth), classified = mean(call != other_class),
                     other_given_class = sum(p$truth == other_class & call != other_class))
    for(cl in real){
      tp = sum(p$truth == cl & call == cl)
      row[[paste0("F1_", cl)]] = if(tp == 0) 0 else 2 * tp / (2 * tp + sum(p$truth != cl & call == cl) + sum(p$truth == cl & call != cl))
    }
    row
  }))
  best = max(scores[[metric]])
  near = scores$k[scores[[metric]] >= best - tolerance]
  recommended = .nearest_grid(stats::median(near), ks)
  if(recommended > smallest){
    warning("the recommended k (", recommended, ") exceeds the smallest class's training samples in a fold (", smallest, ")")
  }
  message("recommended k: ", recommended, " (k within ", tolerance, " of the best ", metric, ": ", paste(near, collapse = ", "), ")")
  plot = ggplot(scores, aes(x = .data$k, y = .data[[metric]])) +
    geom_hline(yintercept = best - tolerance, linetype = "dashed", colour = "grey60") +
    geom_line() + geom_point(aes(colour = .data$k %in% near), size = 2.5) +
    geom_vline(xintercept = recommended, colour = "firebrick") +
    scale_colour_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey40"), guide = "none") +
    labs(x = "k", y = paste("held-out", gsub("_", " ", metric)), title = paste("Recommended k:", recommended)) +
    theme_bw()
  list(recommended_k = recommended, scores = scores, near_best = near, predictions = preds, folds = folds, plot = plot,
       settings = list(ks = ks, metric = metric, tolerance = tolerance, pred_column = pred_column,
                       truth_classes = truth_classes, truth_column = truth_column, other_class = other_class,
                       select_features = select_features, seed = seed))
}
