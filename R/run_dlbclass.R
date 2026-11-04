#' @title Classify DLBCL samples with DLBClass
#'
#' @description Run the DLBClass classifier (Chapuy et al.; the DLBclass-tool
#' neural network models) on one or more genetic sample matrices (GSMs), e.g.
#' from [dlbclass_gsm]. DLBclass-tool is not part of GAMBLR: install it
#' separately (see Details) and point `tool` and `python` at it.
#'
#' @details Setup, once:
#' \preformatted{
#' git clone https://github.com/getzlab/DLBclass-tool.git
#' python3 -m venv DLBclass-tool/venv
#' DLBclass-tool/venv/bin/pip install -r DLBclass-tool/requirements.txt
#' }
#' then give `tool` and `python` (e.g. "DLBclass-tool/venv/bin/python") to
#' each call, or set them for the session with
#' `options(GAMBLR.predict.dlbclass_tool = ..., GAMBLR.predict.dlbclass_python
#' = ...)` or the environment variables `DLBCLASS_TOOL` and `DLBCLASS_PYTHON`.
#' DLBclass-tool pins pandas 1.3.5 and torch 1.10; pandas 2.3 and torch 2.13
#' give the same calls.
#'
#' The GSMs are classified by the bundled script `run_dlbclass.py`, which calls
#' DLBclass-tool's own functions (its feature reduction and classification),
#' so the calls are the tool's. Before running, the GSMs are checked for every
#' feature the classifier uses (read from the tool's `src/format_data.py`).
#'
#' @param gsm A GSM, or a list of them (unnamed ones are named input1,
#' input2, ...): [dlbclass_gsm] output (either
#' `output_format`), or the path of a GSM file in DLBClass's layout (a
#' feature-name column, then one column per sample). CCF, PLOIDY, PURITY and
#' COO rows are ignored, and CN feature names without DLBClass's "X" prefix
#' (e.g. 17P.DEL) get it.
#' @param tool Path of the DLBclass-tool directory. Default: the
#' `GAMBLR.predict.dlbclass_tool` option, else `DLBCLASS_TOOL`.
#' @param python Python with DLBclass-tool's requirements. Default: the
#' `GAMBLR.predict.dlbclass_python` option, else `DLBCLASS_PYTHON`, else
#' "python3".
#' @param min_confidence Calls with Confidence above this get their cluster
#' in the `DLBClass` column; the rest are "Other". Default 0.7.
#'
#' @return For one GSM, a data frame with one row per sample: sample_id, the
#' confidence for each cluster (C1-C5), Confidence (the highest),
#' PredictedCluster and DLBClass. For a list, a list of them with the same
#' names. The
#' DLBclass-tool commit, versions and `min_confidence` are in
#' `attr(result, "dlbclass_run")`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' options(GAMBLR.predict.dlbclass_tool = "~/git/DLBclass-tool",
#'         GAMBLR.predict.dlbclass_python = "~/git/DLBclass-tool/venv/bin/python")
#' gsm = dlbclass_gsm(maf, meta, seg_data = segs)
#' calls = run_dlbclass(gsm)
#' table(calls$DLBClass)
#'
#' # several GSMs (e.g. different copy number inputs), compared
#' calls = run_dlbclass(list(cnv_master = gsm, qc_score = gsm_qc))
#' compare_dlbclass_calls(calls)$summary
#' }
run_dlbclass = function(gsm,
                        tool = NULL,
                        python = NULL,
                        min_confidence = 0.7){
  if(is.null(tool)){
    tool = getOption("GAMBLR.predict.dlbclass_tool", Sys.getenv("DLBCLASS_TOOL"))
  }
  if(is.null(python)){
    python = getOption("GAMBLR.predict.dlbclass_python", Sys.getenv("DLBCLASS_PYTHON"))
    if(identical(python, "")) python = "python3"
  }
  tool = path.expand(tool)
  python = path.expand(python)
  if(identical(tool, "") || !dir.exists(file.path(tool, "src")) ||
     !dir.exists(file.path(tool, "saved_models"))){
    stop("DLBclass-tool not found", if(!identical(tool, "")) paste0(" at ", tool) else "",
         ": clone https://github.com/getzlab/DLBclass-tool and give its directory as tool ",
         "(or set options(GAMBLR.predict.dlbclass_tool = ...)); see ?run_dlbclass")
  }
  script = system.file("python", "run_dlbclass.py", package = "GAMBLR.predict")
  versions = suppressWarnings(system2(python, c(shQuote(script), "--versions"),
                                      stdout = TRUE, stderr = TRUE))
  if(!is.null(attr(versions, "status")) || !any(grepl("^python=", versions))){
    stop("the Python ", python, " can't run DLBclass-tool (it needs pandas and torch; see ?run_dlbclass):\n",
         paste(utils::tail(versions, 3), collapse = "\n"))
  }

  single = is.data.frame(gsm) || is.matrix(gsm) || is.character(gsm)
  gsms = if(single) list(gsm = gsm) else gsm
  # unnamed ones are named input1, input2, ...
  if(is.null(names(gsms))) names(gsms) = rep("", length(gsms))
  unnamed = names(gsms) == ""
  names(gsms)[unnamed] = paste0("input", which(unnamed))
  gsms = lapply(gsms, .dlbclass_gsm_matrix)
  required = .dlbclass_required_features(tool)
  for(n in names(gsms)){
    absent = setdiff(required, colnames(gsms[[n]]))
    if(length(absent) > 0){
      stop(if(single) "the GSM" else paste0("GSM ", n), " lacks ", length(absent),
           " feature(s) the classifier uses: ", paste(utils::head(absent, 20), collapse = ", "),
           if(length(absent) > 20) ", ..." else "")
    }
  }

  dir = tempfile("dlbclass_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  gsm_files = file.path(dir, paste0(seq_along(gsms), ".GSM.tsv"))
  call_files = file.path(dir, paste0(seq_along(gsms), ".classified.tsv"))
  for(i in seq_along(gsms)){
    m = gsms[[i]]
    utils::write.table(data.frame(classifier_name = colnames(m), t(m), check.names = FALSE),
                       gsm_files[i], sep = "\t", quote = FALSE, row.names = FALSE)
  }
  log = suppressWarnings(system2(python, c(shQuote(script), "--tool", shQuote(tool),
                                           as.vector(rbind(shQuote(gsm_files), shQuote(call_files)))),
                                 stdout = TRUE, stderr = TRUE))
  if(!is.null(attr(log, "status")) || !all(file.exists(call_files))){
    stop("DLBclass-tool failed:\n", paste(utils::tail(log, 10), collapse = "\n"))
  }

  commit = suppressWarnings(tryCatch(system2("git", c("-C", shQuote(tool), "rev-parse", "--short", "HEAD"),
                                             stdout = TRUE, stderr = FALSE), error = function(e) NA_character_))
  run = list(tool = tool, commit = if(length(commit) == 1) commit else NA_character_,
             versions = versions[grepl("^python=", versions)], min_confidence = min_confidence)
  calls = lapply(call_files, function(f){
    x = utils::read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
    names(x)[1] = "sample_id"
    x$DLBClass = ifelse(x$Confidence > min_confidence, x$PredictedCluster, "Other")
    attr(x, "dlbclass_run") = run
    x
  })
  names(calls) = names(gsms)
  if(single) calls[[1]] else calls
}

#' @title Compare DLBClass calls
#'
#' @description Compare several sets of DLBClass calls of the same samples
#' (e.g. from GSMs built with different copy number inputs) with the first
#' set, the reference: how often the cluster and the label agree, a
#' confusion table, and the samples whose label changed.
#'
#' @param calls Named list of calls from [run_dlbclass] (or data frames with
#' sample_id, Confidence and PredictedCluster); the first is the reference.
#' @param min_confidence Calls with Confidence above this keep their cluster
#' as their label; the rest are "Other". Default 0.7, as in [run_dlbclass].
#' NULL compares the clusters only.
#'
#' @return A list: `summary` (one row per compared set: n_samples shared,
#' same_cluster (fraction with the same PredictedCluster), same_label (the
#' same label, "Other" included), kappa (Cohen's, of the labels), n_changed,
#' and the fractions of confident calls in the reference and the set),
#' `confusion` (a table per set, reference labels as rows) and `changed` (per
#' set, the samples whose label differs, with both labels and confidences).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' calls = run_dlbclass(list(cnv_master = gsm, qc_score = gsm_qc))
#' cc = compare_dlbclass_calls(calls)
#' cc$summary
#' cc$confusion$qc_score
#' cc$changed$qc_score
#' }
compare_dlbclass_calls = function(calls, min_confidence = 0.7){
  if(!is.list(calls) || is.data.frame(calls) || length(calls) < 2 ||
     is.null(names(calls)) || any(names(calls) == "")){
    stop("give calls as a named list of at least two sets")
  }
  label = function(x){
    cl = as.character(x$PredictedCluster)
    if(!is.null(min_confidence)) cl[!(x$Confidence > min_confidence)] = "Other"
    stats::setNames(cl, x$sample_id)
  }
  ref_name = names(calls)[1]
  ref = calls[[1]]
  ref_label = label(ref)
  per_set = lapply(names(calls)[-1], function(n){
    x = calls[[n]]
    x_label = label(x)
    s = intersect(names(ref_label), names(x_label))
    if(length(s) == 0) stop("no samples in common between ", ref_name, " and ", n)
    ri = match(s, ref$sample_id)
    xi = match(s, x$sample_id)
    same_cluster = as.character(ref$PredictedCluster[ri]) == as.character(x$PredictedCluster[xi])
    changed = data.frame(sample_id = s, ref_label = ref_label[s], label = x_label[s],
                         ref_confidence = ref$Confidence[ri], confidence = x$Confidence[xi],
                         stringsAsFactors = FALSE)
    changed = changed[changed$ref_label != changed$label, , drop = FALSE]
    names(changed)[2:5] = c(ref_name, n, paste0(ref_name, "_confidence"), paste0(n, "_confidence"))
    rownames(changed) = NULL
    levels = sort(unique(c(ref_label[s], x_label[s])))
    list(summary = data.frame(set = n, n_samples = length(s),
                              same_cluster = mean(same_cluster),
                              same_label = mean(ref_label[s] == x_label[s]),
                              kappa = .cohens_kappa(ref_label[s], x_label[s]),
                              n_changed = nrow(changed),
                              confident_reference = mean(ref_label[s] != "Other"),
                              confident = mean(x_label[s] != "Other"),
                              stringsAsFactors = FALSE),
         confusion = table(factor(ref_label[s], levels), factor(x_label[s], levels),
                           dnn = c(ref_name, n)),
         changed = changed)
  })
  names(per_set) = names(calls)[-1]
  summary = do.call(rbind, lapply(per_set, `[[`, "summary"))
  rownames(summary) = NULL
  list(reference = ref_name, min_confidence = min_confidence, summary = summary,
       confusion = lapply(per_set, `[[`, "confusion"),
       changed = lapply(per_set, `[[`, "changed"))
}

# A GSM as a numeric matrix (samples as rows, features as columns) from
# dlbclass_gsm output (either format) or a GSM file; CCF, PLOIDY, PURITY and
# COO rows are dropped and CN features get DLBClass's "X" prefix
.dlbclass_gsm_matrix = function(x){
  if(is.character(x) && length(x) == 1){
    x = utils::read.delim(x, check.names = FALSE, stringsAsFactors = FALSE)
  }
  x = as.data.frame(x, check.names = FALSE)
  if(colnames(x)[1] %in% c("classifier_name", "0", "") || !is.numeric(x[[1]])){
    features = as.character(x[[1]])
    m = t(as.matrix(x[, -1, drop = FALSE]))
    colnames(m) = features
  }else{
    m = as.matrix(x)
  }
  m = m[, !grepl("\\.CCF$", colnames(m)) & !colnames(m) %in% c("PLOIDY", "PURITY", "COO"), drop = FALSE]
  storage.mode(m) = "numeric"
  cn = grepl("^[0-9]+[PQ]", colnames(m))
  colnames(m)[cn] = paste0("X", colnames(m)[cn])
  m
}

# Features the DLBclass-tool classifier reads: the DLBClass features named in
# its src/format_data.py
.dlbclass_required_features = function(tool){
  src = paste(readLines(file.path(tool, "src", "format_data.py"), warn = FALSE), collapse = "\n")
  quoted = gsub('"', "", regmatches(src, gregexpr('"[A-Za-z0-9_.]+"', src))[[1]])
  order = readLines(system.file("extdata", "dlbclass_feature_order.txt", package = "GAMBLR.predict"))
  intersect(order, quoted)
}

# Cohen's kappa of two label vectors (NA when undefined)
.cohens_kappa = function(a, b){
  levels = sort(unique(c(a, b)))
  if(length(levels) < 2) return(NA_real_)
  t = table(factor(a, levels), factor(b, levels))
  n = sum(t)
  observed = sum(diag(t)) / n
  expected = sum(rowSums(t) * colSums(t)) / n^2
  if(expected == 1) return(NA_real_)
  (observed - expected) / (1 - expected)
}
