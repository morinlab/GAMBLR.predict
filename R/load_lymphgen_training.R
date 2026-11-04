#' Load the LymphGen training data as a feature matrix
#'
#' Returns the LymphGen training data (Wright et al. 2020, Cancer Cell;
#' distributed with the LymphGen R code) as a samples x features data frame
#' with metadata, in the layout [make_and_annotate_umap()] expects: one row
#' per sample (row names are sample IDs), one numeric column per feature, and
#' metadata with `sample_id` and `lymphgen` columns. It also returns the
#' features LymphGen selects, with its weights, and every feature it doesn't
#' use.
#'
#' The data aren't distributed with the package: build them from the LymphGen
#' R distribution with `data-raw/lymphgen_training.R` (in the GAMBLR.predict
#' repository), then pass the file as `lymphgen_data` or set
#' `options(GAMBLR.predict.lymphgen_training = <path>)`. A copy placed in
#' `inst/extdata/lymphgen_training.rds` is used when present. The models
#' are those of LymphGen's own `Model9()` run on the training data with every
#' gene available: model 1 (Full) for BN2, EZB, MCD, N1 and ST2, and model 9
#' for A53. The `lymphgen` labels are the training labels (Stage 3 GenClass,
#' with SDT renamed ST2).
#'
#' Some training data are missing: copy number features are NA in the 14
#' samples without copy number, and mutation features are NA in 15 samples for
#' genes their sequencing didn't cover. `complete_only = TRUE` removes samples
#' with any NA in the returned features (29 for the model features).
#' Otherwise, use `na_vals = "to_zero"` in [make_and_annotate_umap()]: its
#' default, `"drop"`, removes every column with an NA.
#'
#' @param feature_set `"model"` (default): one column per feature LymphGen
#'   selected (115 with A53), named subtype:gene:levels, e.g.
#'   `MCD:CDKN2A:HOMDEL|MUTATION`, valued according to `encoding`.
#'   `"atomic"`: one 0/1 column per gene and base feature type (TRUNC,
#'   MUTATION, Synon, HOMDEL, HETLOSS, AMP, GAIN), plus MYD88 L265P and the
#'   BCL2 and BCL6 fusions, present in at least `min_count` samples. The types
#'   are cumulative, as in the training data: MUTATION includes TRUNC, HETLOSS
#'   includes HOMDEL, GAIN includes AMP.
#' @param encoding How model features are valued (ignored for `"atomic"`):
#'   * `"level"`: the model's level: 0 (absent), 1 (core level, e.g. TRUNC or
#'     HOMDEL) or 2+ (wider levels, e.g. MUTATION or HETLOSS). The core level,
#'     LymphGen's strongest evidence, gets the smaller value.
#'   * `"ranked"`: levels reversed so the core level is highest: for a feature
#'     with n levels, level L becomes n - L + 1 (absent = 0). E.g. ID3
#'     (TRUNC, then MUTATION): truncating 2, other coding mutation 1.
#'   * `"presence"`: 1 if present at any level, else 0.
#'   * `"weighted"`: LymphGen's weight for the sample's level minus its weight
#'     when the feature is absent (absent = 0): the feature's contribution to
#'     the sample's log score for that subtype.
#' @param include_a53 Include the A53 model's features (default TRUE).
#' @param min_count For `"atomic"`, the minimum number of training samples
#'   with a feature (default 5).
#' @param complete_only Keep only samples with no missing data in the
#'   returned features (default TRUE).
#' @param sample_links Optional data frame (or path to a TSV) linking training
#'   samples to other IDs, with columns `study_order` (the training sample's
#'   column in the LymphGen training matrix) and `sample_id`; added to the
#'   metadata as `linked_sample_id`.
#' @param lymphgen_data Path to the training data built by
#'   `data-raw/lymphgen_training.R`. Default: the
#'   `GAMBLR.predict.lymphgen_training` option, or else a copy in the
#'   package's `extdata`, if any.
#'
#' @return A list with
#'   * `df`: samples x features data frame; samples are named by their
#'     training column, LG_001 to LG_574.
#'   * `metadata`: `sample_id`, `study_order`, `lymphgen` (training label),
#'     `has_cn` and, with `sample_links`, `linked_sample_id`.
#'   * `feature_weights`: one row per level of each feature LymphGen selected
#'     (`subtype`, `model_feature` = the `df` column name with
#'     `feature_set = "model"`, `gene`, `entrez_id`, `feature_index` = row of
#'     the training annotation, `feature_type`, `level`, `n_levels`), with
#'     `weight` (LymphGen's weight for that level), `absent_weight` (its weight
#'     when the feature is absent), `contribution` (`weight - absent_weight`),
#'     and the number of training samples in (`n_in_subtype`) and outside
#'     (`n_outside`) the subtype with that level.
#'   * `dropped_features`: every other training feature, with the number of
#'     training samples that have it and the `reason` it isn't used:
#'     `"prefilter"` (failed LymphGen's prefilters: too rare, mostly from the
#'     minority data type, plain MYD88, arms or single-copy gains/losses
#'     outside A53, ...) or `"not_selected"` (passed those but not the
#'     per-subtype Fisher test and hierarchy building).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' lg = load_lymphgen_training(encoding = "weighted")
#' u = make_and_annotate_umap(lg$df, lg$metadata)
#' make_umap_scatterplot(u$df)
#' head(lg$feature_weights)
#' }
load_lymphgen_training = function(feature_set = c("model", "atomic"),
                                  encoding = c("level", "ranked", "presence", "weighted"),
                                  include_a53 = TRUE,
                                  min_count = 5,
                                  complete_only = TRUE,
                                  sample_links = NULL,
                                  lymphgen_data = getOption("GAMBLR.predict.lymphgen_training", system.file("extdata", "lymphgen_training.rds",
                                                              package = "GAMBLR.predict"))){
  feature_set = match.arg(feature_set)
  encoding = match.arg(encoding)
  if(!nzchar(lymphgen_data) || !file.exists(lymphgen_data)){
    stop("LymphGen training data not found", if(nzchar(lymphgen_data)) paste0(": ", lymphgen_data) else "", ". ",
         "They aren't distributed with GAMBLR.predict: build them from the LymphGen R distribution with ",
         "data-raw/lymphgen_training.R (in the GAMBLR.predict repository) and pass the file as lymphgen_data, ",
         "or set options(GAMBLR.predict.lymphgen_training = <path>)")
  }
  lg = readRDS(lymphgen_data)
  uber = lg$uber; study = lg$study; models = lg$models[c("Full", "A53")]
  n_feat = lg$fullmat$dim[1]
  ids = sprintf("LG_%03d", study$Study.order)
  label = study$Stage3.genclass
  label[label == "SDT"] = "ST2"
  metadata = data.frame(sample_id = ids, study_order = study$Study.order, lymphgen = label,
                        has_cn = study$CGH == 1, stringsAsFactors = FALSE)
  if(!is.null(sample_links)){
    if(is.character(sample_links)){
      sample_links = utils::read.delim(sample_links, stringsAsFactors = FALSE)
    }
    if(!all(c("study_order", "sample_id") %in% colnames(sample_links))){
      stop("sample_links needs study_order and sample_id columns")
    }
    metadata$linked_sample_id = sample_links$sample_id[match(metadata$study_order, sample_links$study_order)]
  }

  # rows of the training matrix (0/1/NA), rebuilt from the stored cell indices
  feature_rows = function(rows){
    out = matrix(0, length(rows), length(ids))
    pos = match(seq_len(n_feat), rows)
    for(what in c("ones", "na")){
      cell = lg$fullmat[[what]]
      r = pos[(cell - 1) %% n_feat + 1]
      hit = !is.na(r)
      out[cbind(r[hit], (cell[hit] - 1) %/% n_feat + 1)] = if(what == "ones") 1 else NA
    }
    out
  }

  full_types = sort(unique(label[!label %in% c("Other", "A53")]))
  clean_gene = function(gene){
    gene[gene == "BCL6FusionCatcherORInHouse"] = "BCL6"
    gene[gene == "BCL2TransFISHMUTFusion"] = "BCL2"
    gene
  }
  # model feature rows (dropping the empty padding rows) and, per row, the
  # training annotation rows of its levels: out (core), ot1, ot1.1 (wider)
  level_index = function(model){
    mt = model$Mod.type
    keep = !is.na(mt$out) & mt$out > 0
    list(keep = keep, idx = as.matrix(mt[keep, intersect(c("out", "ot1", "ot1.1"), colnames(mt)), drop = FALSE]))
  }
  model_samples = list(rep(TRUE, length(ids)), study$CGH == 1)  # A53 is trained on samples with copy number
  model_subtypes = list(full_types, "A53")
  use_models = if(include_a53) 1:2 else 1

  feature_weights = NULL
  for(k in use_models){
    model = models[[k]]
    li = level_index(model); mt = model$Mod.type[li$keep, ]
    w = as.matrix(model$Mod.scale[li$keep, , drop = FALSE])
    sd = as.matrix(model$startdat[li$keep, , drop = FALSE])
    lab_k = label[model_samples[[k]]]
    for(r in seq_len(nrow(mt))){
      fi = li$idx[r, ]; fi = fi[!is.na(fi) & fi > 0]
      st = model_subtypes[[k]][mt$i[r]]
      inside = lab_k == st
      feature_weights = rbind(feature_weights, data.frame(
        subtype = st,
        model_feature = paste(st, clean_gene(uber$Gene.Symbol[fi[1]]),
                              paste(uber$Feature.type[fi], collapse = "|"), sep = ":"),
        gene = clean_gene(uber$Gene.Symbol[fi]), entrez_id = uber$ENTREZ.ID[fi],
        feature_index = fi, feature_type = uber$Feature.type[fi],
        level = seq_along(fi), n_levels = length(fi),
        weight = w[r, seq_along(fi) + 1], absent_weight = w[r, 1],
        contribution = w[r, seq_along(fi) + 1] - w[r, 1],
        n_in_subtype = vapply(seq_along(fi), function(L) sum(sd[r, inside] == L, na.rm = TRUE), numeric(1)),
        n_outside = vapply(seq_along(fi), function(L) sum(sd[r, !inside] == L, na.rm = TRUE), numeric(1)),
        stringsAsFactors = FALSE))
    }
  }
  rownames(feature_weights) = NULL

  n_samples = tabulate((lg$fullmat$ones - 1) %% n_feat + 1, nbins = n_feat)
  dropped = !uber$Feature.Index %in% feature_weights$feature_index
  dropped_features = data.frame(feature_index = uber$Feature.Index[dropped],
                                gene = uber$Gene.Symbol[dropped], entrez_id = uber$ENTREZ.ID[dropped],
                                feature_type = uber$Feature.type[dropped],
                                n_samples = n_samples[dropped],
                                reason = ifelse(lg$prefilter[dropped], "not_selected", "prefilter"),
                                stringsAsFactors = FALSE)

  if(feature_set == "atomic"){
    base = c("TRUNC", "MUTATION", "Synon", "HOMDEL", "HETLOSS", "AMP", "GAIN")
    keep = uber$Feature.type %in% base & uber$ENTREZ.ID != 1000000004
    keep = keep | (uber$ENTREZ.ID == 1000000004 & uber$Feature.type == "MUTATION")  # MYD88 L265P
    keep = keep | uber$Feature.Index %in% c(3534, 3619)                                 # BCL2 and BCL6 fusions
    keep = keep & n_samples >= min_count
    rows = which(keep)
    m = feature_rows(rows)
    nm = paste(uber$Gene.Symbol[rows], uber$Feature.type[rows], sep = ":")
    nm[uber$ENTREZ.ID[rows] == 1000000004] = "MYD88:L265P"
    nm[uber$Feature.Index[rows] == 3534] = "BCL2:FUSION"
    nm[uber$Feature.Index[rows] == 3619] = "BCL6:FUSION"
    rownames(m) = nm
  }else{
    m = NULL
    for(k in use_models){
      model = models[[k]]
      li = level_index(model); keep = li$keep; idx = li$idx
      lev = apply(idx, 1, function(r) paste(uber$Feature.type[r[!is.na(r) & r > 0]], collapse = "|"))
      gene = clean_gene(uber$Gene.Symbol[model$Mod.type$out[keep]])
      nm = paste(model_subtypes[[k]][model$Mod.type$i[keep]], gene, lev, sep = ":")
      lv = as.matrix(model$startdat[keep, , drop = FALSE])
      if(encoding == "ranked"){
        n_levels = rowSums(!is.na(idx) & idx > 0)
        lv = ifelse(lv > 0, pmax(n_levels - lv + 1, 1), lv)
      }else if(encoding == "presence"){
        lv = (lv > 0) + 0
      }else if(encoding == "weighted"){
        # Mod.scale column L + 1 is the weight for level L (column 1: absent)
        w = as.matrix(model$Mod.scale[keep, , drop = FALSE])
        lv = matrix(w[cbind(rep(seq_len(nrow(lv)), ncol(lv)), as.vector(lv) + 1)] - w[, 1],
                    nrow(lv), dimnames = dimnames(lv))
      }
      d = matrix(NA_real_, sum(keep), length(ids), dimnames = list(nm, ids))
      d[, model_samples[[k]]] = lv
      m = rbind(m, d)
    }
  }

  df = as.data.frame(t(m), check.names = FALSE)
  rownames(df) = ids
  if(complete_only){
    complete = rowSums(is.na(df)) == 0
    message(sprintf("keeping %d of %d samples with complete data for these features", sum(complete), nrow(df)))
    df = df[complete, , drop = FALSE]
    metadata = metadata[complete, ]
  }
  list(df = df, metadata = metadata, feature_weights = feature_weights,
       dropped_features = dropped_features)
}
