#' @title Compare two sets of LymphGen inputs
#'
#' @description Find what differs between two [lymphgen_input] results for
#' the same samples, e.g. built from different `seg_data` or `cnv_args`, and
#' summarise it layer by layer: sample annotation, gene list, mutations, gene
#' and arm copy number, and the LymphGen features those produce, by subgroup.
#' With the LymphGen calls of both ([run_lymphgen]), the differing features of
#' each sample whose call changed are listed, to show what drove the change.
#'
#' @details Copy number states are compared per sample and gene (or arm) over
#' the samples with copy number in both inputs: a change is "a only" or "b
#' only" (called in one input only), "level" (HETLOSS vs HOMDEL, or GAIN vs
#' AMP) or "direction" (a loss in one and a gain in the other). Features are
#' compared in the inputs' feature matrices (1 = present); a feature
#' assessed in one input but not the other (NA, e.g. copy number features
#' without copy number) is counted as "unassessed". Each feature's subgroup
#' and data type come from the LymphGen features of the GAMBLR.utils
#' classifier registry ([GAMBLR.utils::get_classifier_features]).
#'
#' @param a,b Two [lymphgen_input] results.
#' @param names Names for `a` and `b` in the output. Default c("a", "b").
#' @param calls Optional LymphGen calls for `a` and `b`: a list of two
#' [run_lymphgen] results, in the same order.
#'
#' @return A list (printed as a summary):
#' \itemize{
#'   \item `overview`: one row per layer with the number of differences.
#'   \item `samples`: per shared sample, its differing mutations, copy number
#'     states and features, the annotation differences and, with `calls`, both
#'     calls.
#'   \item `annotation`: samples whose Copy.Number or translocation status
#'     differs.
#'   \item `gene_list`: Entrez IDs only in `a` or only in `b`.
#'   \item `mutations`: mutation rows only in one input (sample, gene, Type,
#'     Location).
#'   \item `copy_number`: differing gene and arm copy number states, with the
#'     kind of difference.
#'   \item `features`: differing feature calls (sample, feature, subgroup,
#'     data type, both values, kind).
#'   \item `by_feature` and `by_subgroup`: those counted per feature and per
#'     subgroup and data type.
#'   \item `changed_calls`: with `calls`, each sample whose call changed, with
#'     both calls and its differing features in the subgroups involved.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' default = lymphgen_input(maf, meta, seg_data = segs)
#' strict = lymphgen_input(maf, meta, seg_data = segs,
#'                         cnv_args = list(single_thresh = 0.7))
#' cmp = compare_lymphgen_inputs(default, strict, names = c("default", "strict"),
#'                               calls = list(run_lymphgen(default), run_lymphgen(strict)))
#' cmp
#' cmp$by_subgroup
#' cmp$changed_calls
#' }
compare_lymphgen_inputs = function(a, b, names = c("a", "b"), calls = NULL){
  for(x in list(a, b)){
    if(!is.list(x) || !all(c("mutation_flat", "gene_list", "sample_annotation") %in% base::names(x))){
      stop("a and b must be lymphgen_input() results")
    }
  }
  registry = GAMBLR.utils::get_classifier_features("LymphGen")
  # feature matrix columns are LymphGen's feature names (classifier_feature);
  # a feature can belong to several subgroups or data types
  collapse = function(x) paste(sort(unique(stats::na.omit(x))), collapse = "/")
  subgroup_of = tapply(registry$group, registry$classifier_feature, collapse)
  data_type_of = tapply(registry$data_type, registry$classifier_feature, collapse)
  symbols = stats::setNames(registry$gene, registry$entrez)
  symbols = symbols[!is.na(base::names(symbols)) & !duplicated(base::names(symbols))]
  gene_name = function(id) ifelse(is.na(symbols[as.character(id)]), as.character(id), symbols[as.character(id)])

  annot_a = as.data.frame(a$sample_annotation)
  annot_b = as.data.frame(b$sample_annotation)
  samples = intersect(annot_a$Sample.ID, annot_b$Sample.ID)
  if(length(samples) == 0) stop("a and b have no samples in common")
  only_samples = list(setdiff(annot_a$Sample.ID, samples), setdiff(annot_b$Sample.ID, samples))
  if(length(unlist(only_samples)) > 0){
    message("comparing the ", length(samples), " shared samples (", length(only_samples[[1]]), " only in ",
            names[1], ", ", length(only_samples[[2]]), " only in ", names[2], ")")
  }

  # sample annotation: copy number availability and translocations
  annot_a = annot_a[match(samples, annot_a$Sample.ID), , drop = FALSE]
  annot_b = annot_b[match(samples, annot_b$Sample.ID), , drop = FALSE]
  annot_cols = intersect(c("Copy.Number", "BCL2.transloc", "BCL6.transloc"),
                         intersect(colnames(annot_a), colnames(annot_b)))
  annotation = do.call(rbind, lapply(annot_cols, function(col){
    va = as.character(annot_a[[col]]); vb = as.character(annot_b[[col]])
    # row by row: one NA and not the other, or both known and different
    differ = (is.na(va) != is.na(vb)) | (!is.na(va) & !is.na(vb) & va != vb)
    if(!any(differ)) return(NULL)
    data.frame(sample_id = samples[differ], field = col, a = va[differ], b = vb[differ],
               stringsAsFactors = FALSE)
  }))
  if(is.null(annotation)) annotation = data.frame(sample_id = character(0), field = character(0),
                                                  a = character(0), b = character(0))

  # gene list
  genes_a = as.character(a$gene_list[[1]]); genes_b = as.character(b$gene_list[[1]])
  gene_list = list(only_a = setdiff(genes_a, genes_b), only_b = setdiff(genes_b, genes_a))

  # mutation rows only in one input
  mut_key = function(m){
    m = as.data.frame(m)
    m = m[m$Sample %in% samples, , drop = FALSE]
    if(!"Location" %in% colnames(m)) m$Location = -1
    m$key = paste(m$Sample, m$ENTREZ.ID, m$Type, m$Location)
    m
  }
  ma = mut_key(a$mutation_flat); mb = mut_key(b$mutation_flat)
  only_in = function(x, y, name){
    x = x[!x$key %in% y$key, c("Sample", "ENTREZ.ID", "Type", "Location"), drop = FALSE]
    x$in_input = rep(name, nrow(x))
    x
  }
  mutations = rbind(only_in(ma, mb, names[1]), only_in(mb, ma, names[2]))
  mutations = data.frame(sample_id = mutations$Sample, gene = gene_name(mutations$ENTREZ.ID),
                         mutations[, c("ENTREZ.ID", "Type", "Location", "in_input")],
                         stringsAsFactors = FALSE)

  # gene and arm copy number states, over samples with copy number in both
  with_cn = function(annot) annot$Sample.ID[as.character(annot$Copy.Number) %in% c("1", "TRUE")]
  cn_samples = intersect(with_cn(annot_a), with_cn(annot_b))
  cn_states = function(x, region_col){
    if(is.null(x) || nrow(x) == 0) return(data.frame(Sample = character(0), region = character(0), Type = character(0)))
    x = as.data.frame(x)
    x = x[x$Sample %in% cn_samples, , drop = FALSE]
    data.frame(Sample = x$Sample, region = as.character(x[[region_col]]), Type = as.character(x$Type),
               stringsAsFactors = FALSE)
  }
  direction = function(type) ifelse(type %in% c("HETLOSS", "HOMDEL"), "loss",
                                    ifelse(type %in% c("GAIN", "AMP"), "gain", NA))
  compare_cn = function(sa, sb, level){
    both = merge(sa, sb, by = c("Sample", "region"), all = TRUE, suffixes = c("_a", "_b"))
    same = !is.na(both$Type_a) & !is.na(both$Type_b) & both$Type_a == both$Type_b
    both = both[!same, , drop = FALSE]
    if(nrow(both) == 0) return(NULL)
    kind = ifelse(is.na(both$Type_b), paste(names[1], "only"),
                  ifelse(is.na(both$Type_a), paste(names[2], "only"),
                         ifelse(direction(both$Type_a) == direction(both$Type_b), "level", "direction")))
    data.frame(sample_id = both$Sample, level = level,
               region = if(level == "gene") gene_name(both$region) else both$region,
               a = both$Type_a, b = both$Type_b, kind = kind, stringsAsFactors = FALSE)
  }
  copy_number = rbind(compare_cn(cn_states(a$copy_flat, "ENTREZ.ID"), cn_states(b$copy_flat, "ENTREZ.ID"), "gene"),
                      compare_cn(cn_states(a$arm, "Arm"), cn_states(b$arm, "Arm"), "arm"))
  if(is.null(copy_number)) copy_number = data.frame(sample_id = character(0), level = character(0),
                                                    region = character(0), a = character(0),
                                                    b = character(0), kind = character(0))

  # features
  features = NULL
  if(!is.null(a$feature_matrix) && !is.null(b$feature_matrix)){
    fa = as.matrix(a$feature_matrix); fb = as.matrix(b$feature_matrix)
    fs = intersect(rownames(fa), samples); fs = intersect(fs, rownames(fb))
    ff = intersect(colnames(fa), colnames(fb))
    A = fa[fs, ff, drop = FALSE]; B = fb[fs, ff, drop = FALSE]
    differ = (is.na(A) != is.na(B)) | (!is.na(A) & !is.na(B) & A != B)
    at = which(differ, arr.ind = TRUE)
    if(nrow(at) > 0){
      feat = ff[at[, 2]]
      va = A[at]; vb = B[at]
      features = data.frame(sample_id = fs[at[, 1]], feature = feat,
                            subgroup = as.character(subgroup_of[feat]),
                            data_type = as.character(data_type_of[feat]),
                            a = va, b = vb,
                            kind = ifelse(is.na(va) | is.na(vb), "unassessed",
                                          ifelse(va > 0 & !vb > 0, paste(names[1], "only"),
                                                 ifelse(vb > 0 & !va > 0, paste(names[2], "only"), "level"))),
                            stringsAsFactors = FALSE)
    }
  }
  if(is.null(features)) features = data.frame(sample_id = character(0), feature = character(0),
                                              subgroup = character(0), data_type = character(0),
                                              a = numeric(0), b = numeric(0), kind = character(0))
  by_feature = if(nrow(features) > 0){
    f = stats::aggregate(list(n = rep(1, nrow(features))),
                         features[, c("feature", "subgroup", "data_type", "kind")], sum)
    f = stats::reshape(f, idvar = c("feature", "subgroup", "data_type"), timevar = "kind", direction = "wide")
    colnames(f) = sub("^n\\.", "", colnames(f))
    f[is.na(f)] = 0
    f$total = rowSums(f[, setdiff(colnames(f), c("feature", "subgroup", "data_type")), drop = FALSE])
    f = f[order(-f$total), , drop = FALSE]
    rownames(f) = NULL
    f
  } else NULL
  by_subgroup = if(nrow(features) > 0){
    key = paste(features$subgroup, features$data_type, sep = "\r")
    first = !duplicated(key)
    s = data.frame(subgroup = features$subgroup[first], data_type = features$data_type[first],
                   n_differences = as.vector(table(key)[key[first]]),
                   n_samples = as.vector(tapply(features$sample_id, key, function(x) length(unique(x)))[key[first]]),
                   stringsAsFactors = FALSE)
    s = s[order(-s$n_differences), , drop = FALSE]
    rownames(s) = NULL
    s
  } else NULL

  # per sample
  count_by = function(x) as.vector(table(factor(x$sample_id, levels = samples)))
  per_sample = data.frame(sample_id = samples,
                          annotation = count_by(annotation),
                          mutations = count_by(mutations),
                          copy_number = count_by(copy_number),
                          features = count_by(features),
                          stringsAsFactors = FALSE)

  # calls: which differing features sit in the subgroups of the changed calls
  changed_calls = NULL
  if(!is.null(calls)){
    if(length(calls) != 2) stop("calls must be a list of two run_lymphgen() results (for a and b)")
    call_of = function(x) stats::setNames(as.character(x$LymphGen), x$sample_id)
    ca = call_of(calls[[1]]); cb = call_of(calls[[2]])
    per_sample[[paste0("call_", names[1])]] = ca[samples]
    per_sample[[paste0("call_", names[2])]] = cb[samples]
    moved = samples[!is.na(ca[samples]) & !is.na(cb[samples]) & ca[samples] != cb[samples]]
    changed_calls = do.call(rbind, lapply(moved, function(s){
      groups = setdiff(unique(unlist(strsplit(c(ca[s], cb[s]), "/"))), "Other")
      in_groups = vapply(strsplit(features$subgroup, "/"), function(g) any(g %in% groups), TRUE)
      f = features[features$sample_id == s & in_groups, , drop = FALSE]
      data.frame(sample_id = s, a = ca[s], b = cb[s],
                 n_features = nrow(f),
                 features = if(nrow(f) > 0) paste0(f$feature, " (", f$kind, ")", collapse = "; ") else "",
                 stringsAsFactors = FALSE)
    }))
    if(!is.null(changed_calls)){
      colnames(changed_calls)[2:3] = names
      rownames(changed_calls) = NULL
    }
  }
  per_sample = per_sample[order(-(per_sample$features + per_sample$copy_number + per_sample$mutations)), ,
                          drop = FALSE]
  rownames(per_sample) = NULL

  # features first (what LymphGen scores), then the input layers behind them
  overview = data.frame(
    layer = c("features", "annotation", "gene_list", "mutations", "copy_number (gene)", "copy_number (arm)"),
    n_differences = c(nrow(features), nrow(annotation), length(unlist(gene_list)), nrow(mutations),
                      sum(copy_number$level == "gene"), sum(copy_number$level == "arm")),
    n_samples = c(length(unique(features$sample_id)), length(unique(annotation$sample_id)), NA,
                  length(unique(mutations$sample_id)),
                  length(unique(copy_number$sample_id[copy_number$level == "gene"])),
                  length(unique(copy_number$sample_id[copy_number$level == "arm"]))),
    stringsAsFactors = FALSE)

  out = list(names = names, n_samples = length(samples), n_cn_samples = length(cn_samples),
             overview = overview, samples = per_sample, annotation = annotation, gene_list = gene_list,
             mutations = mutations, copy_number = copy_number, features = features,
             by_feature = by_feature, by_subgroup = by_subgroup, changed_calls = changed_calls)
  class(out) = "lymphgen_input_comparison"
  out
}

#' @export
print.lymphgen_input_comparison = function(x, ...){
  cat(sprintf("LymphGen inputs: %s vs %s, %d shared samples (%d with copy number in both)\n\n",
              x$names[1], x$names[2], x$n_samples, x$n_cn_samples))
  print(x$overview, row.names = FALSE)
  if(!is.null(x$by_subgroup)){
    cat("\nFeature differences by subgroup:\n")
    print(x$by_subgroup, row.names = FALSE)
    cat("\nMost often differing features:\n")
    print(utils::head(x$by_feature, 10), row.names = FALSE)
  }
  if(nrow(x$copy_number) > 0){
    cat("\nCopy number input differences by kind:\n")
    print(table(level = x$copy_number$level, kind = x$copy_number$kind))
  }
  if(!is.null(x$changed_calls)){
    cat(sprintf("\n%d sample(s) changed call; see $changed_calls\n", nrow(x$changed_calls)))
    print(utils::head(x$changed_calls[, 1:4], 10), row.names = FALSE)
  }
  invisible(x)
}
