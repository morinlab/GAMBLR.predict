#' @title Build the LymphGen input files
#'
#' @description Build the inputs for LymphGen 2 (Run_Lymphgen.R: the mutation
#' flat file, gene list, copy number flat file, arm file and sample
#' annotation) from GAMBLR pipeline outputs, using the LymphGen features of
#' the GAMBLR.utils classifier registry
#' ([GAMBLR.utils::get_classifier_features]) for the Synon windows, copy
#' number features and translocations.
#'
#' @details
#' \itemize{
#'   \item Mutation flat file (Sample, ENTREZ.ID, Type, Location): one row
#'     per mutation, with its grch37 start as Location, so LymphGen applies
#'     its own position rules (NOTCH1/NOTCH2 PEST, CD79B and the Synon
#'     windows). Types: TRUNC (truncating classes, as in
#'     [GAMBLR.utils::classify_mutations]), L265P (MYD88 p.L265P), MUTATION
#'     (missense, in-frame and nonstop mutations) and Synon (silent and non-coding mutations,
#'     e.g. 5' UTR, intron, splice region and flanking, in a Synon window).
#'     LymphGen counts TRUNC as MUTATION and both as Synon, so each mutation
#'     has one row with its most specific type.
#'   \item Gene list: the Entrez IDs of the genes sequenced. LymphGen refits
#'     its model to these genes, so they must include the unmutated genes.
#'   \item Copy number flat file (Sample, ENTREZ.ID, Type) and arm file
#'     (Sample, Arm, Type): the LymphGen gene and arm copy number features
#'     of the registry (gene features count changes of 30 Mb or less), as
#'     HETLOSS or HOMDEL and GAIN or AMP; LymphGen counts HOMDEL as HETLOSS
#'     and AMP as GAIN. Only for samples with copy number. The levels follow
#'     LymphGen's copy number categories: a change of at least half a copy
#'     is HETLOSS or GAIN, and of at least 1.5 copies HOMDEL or AMP (so
#'     HOMDEL is copy number 0.5 or less and AMP more than 3.5, from a
#'     baseline of 2). As in LymphGen's array training data, changes are
#'     relative to the sample: log ratios centred on the sample's median
#'     and converted to copies of a diploid genome, so in a tetraploid
#'     tumour two copies is a HETLOSS rather than a HOMDEL.
#'   \item Feature matrix: the same mutations, copy number and
#'     translocations called as the registry's LymphGen features, applying
#'     the position rules LymphGen applies (see
#'     [GAMBLR.utils::assemble_classifier_features]); for use in R, e.g. for
#'     DLBCLone. Use [GAMBLR.utils::get_classifier_features] with
#'     `submodel` to select one submodel's features.
#'   \item Sample annotation (Sample.ID, Copy.Number, BCL2.transloc,
#'     BCL6.transloc): Copy.Number is 1 for samples with copy number;
#'     translocations are 1, 0, or NA when unknown, from `sv_data` and/or
#'     the metadata columns in `sv_from_metadata`.
#' }
#'
#' @param maf_data MAF data frame (grch37) of the samples' mutations,
#' including silent ones, e.g. from `get_coding_ssm(..., include_silent =
#' TRUE)`. Non-coding mutations can be in here or in `noncoding_maf`.
#' @param these_samples_metadata Metadata with one row per sample_id; defines
#' the samples (mRNA-seq samples are dropped, with a warning). Optional columns `cn_available` (logical; otherwise samples
#' in `seg_data` have copy number) and those named in `sv_from_metadata`.
#' @param noncoding_maf Optional MAF data frame of non-coding mutations
#' (e.g. 5' UTR and intronic, from `get_ssm_by_regions`) for the Synon
#' feature. Without non-coding mutations, Synon counts silent mutations
#' only.
#' @param seg_data Optional segmented copy number data (columns ID, chrom,
#' start, end, log.ratio) for the copy number features.
#' @param sv_data Optional structural variants (as for
#' [GAMBLR.utils::assemble_classifier_features]) for the BCL2 and BCL6
#' translocations.
#' @param sv_from_metadata Named vector mapping BCL2 and BCL6 to metadata
#' columns with their translocation status. Default `c(BCL2 = "bcl2_ba",
#' BCL6 = "bcl6_ba")`; NULL to use `sv_data` only.
#' @param sv_positive_values Metadata values that indicate a translocation.
#' Default "POS".
#' @param genes Genes sequenced (symbols or Entrez IDs) for the gene list,
#' e.g. a panel's genes. Default NULL lists every gene, as for genomes and
#' exomes.
#' @param synon_regions Where silent and non-coding mutations count as
#' Synon: "registry" (default), the Synon windows of the registry's LymphGen
#' features, which are those LymphGen 2 uses with every gene sequenced; or
#' "gene_tss", the first 4 kb from any active TSS of every gene in the gene
#' list (see [GAMBLR.utils::subset_ssm_by_tss]), for panels, where LymphGen
#' can select other Synon features, and then applies its own windows.
#' @param protein_change_column Column of the MAFs with the protein change.
#' Default NULL uses Protein_Change, or HGVSp_Short if there is none.
#' @param cnv_args Optional named list of further arguments for
#' [GAMBLR.utils::segmented_data_to_cnv_status], e.g. `blacklist`, or to
#' override the LymphGen defaults (`value_scale`, `centering`,
#' `single_thresh = 0.5`, `high_thresh = 1.5`).
#' @param output_dir Optional directory to write the files to, as
#' `<prefix>_mutation_flat.tsv`, `_gene_list.txt`, `_copy_flat.tsv`,
#' `_arm.tsv` and `_sample_annotation.tsv` (copy number files only with
#' `seg_data`).
#' @param prefix File name prefix for `output_dir`. Default "lymphgen".
#'
#' @return A list of data frames: mutation_flat, gene_list, copy_flat and
#' arm (NULL without `seg_data`), sample_annotation, and feature_matrix
#' (samples as rows, the LymphGen features of every submodel as columns, as
#' from [GAMBLR.utils::classifier_feature_matrix]: 1 or 0, NA where not
#' assessed, e.g. copy number features without `seg_data`; not written to
#' `output_dir`). With `output_dir`, returned invisibly, with the file paths
#' as attribute "files".
#'
#' @export
#'
#' @examples
#' \dontrun{
#' meta = get_gambl_metadata() %>% dplyr::filter(pathology == "DLBCL")
#' maf = get_coding_ssm(these_samples_metadata = meta, include_silent = TRUE)
#' synon_regions = get_classifier_features("LymphGen") %>%
#'   dplyr::filter(position_rule == "tss4kb") %>%
#'   dplyr::distinct(chrom, start, end)
#' noncoding = get_ssm_by_regions(regions_bed = synon_regions,
#'                                these_samples_metadata = meta,
#'                                streamlined = FALSE)
#' segs = get_cn_segments(these = meta, projection = "grch37")
#' lymphgen_input(maf, meta, noncoding_maf = noncoding, seg_data = segs,
#'                output_dir = "lymphgen_run", prefix = "dlbcl")
#' # then: Rscript Run_Lymphgen.R -m lymphgen_run/dlbcl_mutation_flat.tsv \
#' #   -g lymphgen_run/dlbcl_gene_list.txt -s lymphgen_run/dlbcl_sample_annotation.tsv \
#' #   -c lymphgen_run/dlbcl_copy_flat.tsv -a lymphgen_run/dlbcl_arm.tsv -o result.tsv
#' }
lymphgen_input = function(maf_data,
                          these_samples_metadata,
                          noncoding_maf = NULL,
                          seg_data = NULL,
                          sv_data = NULL,
                          sv_from_metadata = c(BCL2 = "bcl2_ba", BCL6 = "bcl6_ba"),
                          sv_positive_values = "POS",
                          genes = NULL,
                          synon_regions = "registry",
                          protein_change_column = NULL,
                          cnv_args = list(),
                          output_dir = NULL,
                          prefix = "lymphgen"){
  if(!synon_regions %in% c("registry", "gene_tss")){
    stop("synon_regions must be \"registry\" or \"gene_tss\"")
  }
  if(!"sample_id" %in% colnames(these_samples_metadata)){
    stop("these_samples_metadata needs a sample_id column")
  }
  these_samples_metadata = GAMBLR.utils::drop_mrna_samples(
    as.data.frame(these_samples_metadata), "lymphgen_input")
  samples = unique(as.character(these_samples_metadata$sample_id))
  registry = GAMBLR.utils::get_classifier_features("LymphGen")

  # Entrez IDs: the registry's (LymphGen's own) first, then gene_tss
  tss = GAMBLR.utils::gene_tss
  tss = tss[tss$genome_build == "grch37" & tss$primary, , drop = FALSE]
  known = unique(rbind(registry[!is.na(registry$gene) & !is.na(registry$entrez),
                                c("gene", "entrez")],
                       tss[, c("gene", "entrez")]))
  known = known[!duplicated(known$gene), , drop = FALSE]
  to_entrez = function(symbol) as.integer(known$entrez[match(symbol, known$gene)])

  # gene list
  if(is.null(genes)){
    gene_ids = unique(c(as.integer(tss$entrez), to_entrez(registry$gene)))
  }else{
    genes = as.character(genes)
    is_id = grepl("^[0-9]+$", genes)
    gene_ids = c(as.integer(genes[is_id]), to_entrez(genes[!is_id]))
    unknown = genes[!is_id][is.na(to_entrez(genes[!is_id]))]
    if(length(unknown) > 0){
      message("no Entrez ID for ", length(unknown), " genes, left out of the gene list: ",
              paste(utils::head(unknown, 10), collapse = ", "))
    }
  }
  gene_ids = sort(unique(gene_ids[!is.na(gene_ids)]))

  # mutations -------------------------------------------------------------
  maf = as.data.frame(maf_data)
  if(!is.null(noncoding_maf)){
    noncoding_maf = as.data.frame(noncoding_maf)
    shared = intersect(colnames(maf), colnames(noncoding_maf))
    maf = rbind(maf[, shared, drop = FALSE], noncoding_maf[, shared, drop = FALSE])
  }
  maf = maf[maf$Tumor_Sample_Barcode %in% samples, , drop = FALSE]
  key_cols = intersect(c("Tumor_Sample_Barcode", "Chromosome", "Start_Position",
                         "Reference_Allele", "Tumor_Seq_Allele2"), colnames(maf))
  maf = maf[!duplicated(maf[, key_cols]), , drop = FALSE]
  maf = GAMBLR.utils::classify_mutations(maf, protein_change_column)
  # every sample's mutations, for the feature matrix
  all_maf = maf
  maf$entrez = to_entrez(maf$Hugo_Symbol)
  maf = maf[!is.na(maf$entrez) & maf$entrez %in% gene_ids, , drop = FALSE]
  # coding classes as in the registry's LymphGen MUTATION features; everything
  # else (silent, UTR, intron, splice region, flanking, including 5' UTR start
  # gains) is non-coding and can count as Synon
  coding = maf$mutation_class %in% c("trunc", "missense", "inframe", "nonstop")
  l265p = maf$Hugo_Symbol == "MYD88" &
    (maf$protein_change %in% "p.L265P" | maf$Start_Position == 38182641)
  type = ifelse(l265p, "L265P",
                ifelse(coding & maf$mutation_class %in% "trunc", "TRUNC",
                       ifelse(coding, "MUTATION", "Synon")))

  # Synon: silent and non-coding mutations in a Synon window
  synon = which(type == "Synon")
  in_window = rep(FALSE, length(synon))
  maf_chrom = sub("^chr", "", as.character(maf$Chromosome))
  if(synon_regions == "registry"){
    windows = unique(registry[registry$position_rule %in% "tss4kb",
                              c("gene", "chrom", "start", "end")])
    for(i in seq_len(nrow(windows))){
      in_window = in_window | (maf$Hugo_Symbol[synon] == windows$gene[i] &
                                 maf_chrom[synon] == windows$chrom[i] &
                                 maf$Start_Position[synon] >= windows$start[i] &
                                 maf$Start_Position[synon] <= windows$end[i])
    }
  }else if(length(synon) > 0){
    candidates = maf[synon, , drop = FALSE]
    candidates$row = synon
    near = GAMBLR.utils::subset_ssm_by_tss(candidates,
                                           genes = unique(candidates$Hugo_Symbol),
                                           genome_build = "grch37",
                                           noncoding_only = FALSE)
    in_window = synon %in% near$row[near$tss_gene == near$Hugo_Symbol]
  }
  keep = type != "Synon"
  keep[synon[in_window]] = TRUE
  mutation_flat = data.frame(Sample = maf$Tumor_Sample_Barcode[keep],
                             ENTREZ.ID = maf$entrez[keep],
                             Type = type[keep],
                             Location = maf$Start_Position[keep],
                             stringsAsFactors = FALSE)
  mutation_flat = mutation_flat[order(match(mutation_flat$Sample, samples),
                                      mutation_flat$ENTREZ.ID,
                                      mutation_flat$Location), , drop = FALSE]
  rownames(mutation_flat) = NULL

  # registry features: copy number, translocations and the feature matrix --
  cn_features = registry[registry$data_type == "cnv", , drop = FALSE]
  sv_features = registry[registry$data_type == "sv" &
                           registry$gene %in% c("BCL2", "BCL6"), , drop = FALSE]
  # LymphGen's copy number categories (see details)
  cnv_args = utils::modifyList(list(value_scale = "log_ratio", centering = "median",
                                    single_thresh = 0.5, high_thresh = 1.5),
                               cnv_args)
  wanted = registry[!duplicated(registry$feature), , drop = FALSE]
  # non-coding Synon features are assessed only when there are non-coding
  # mutations (e.g. not from a coding-only MAF)
  has_noncoding = !is.null(noncoding_maf) ||
    any(all_maf$mutation_class %in% c("utr5", "noncoding"))
  feature_table = suppressMessages(GAMBLR.utils::assemble_classifier_features(
    these_samples_metadata,
    features = wanted[, setdiff(colnames(wanted),
                                c("classifier", "classifier_feature", "group",
                                  "feature_order", "value", "min_level",
                                  "submodels")), drop = FALSE],
    maf_data = all_maf,
    noncoding_maf = if(has_noncoding) all_maf,
    seg_data = seg_data,
    sv_data = sv_data,
    sv_from_metadata = sv_from_metadata,
    sv_positive_values = sv_positive_values,
    genome_build = "grch37",
    cnv_args = cnv_args
  ))
  status = function(feature){
    calls = feature_table[feature_table$feature == feature, , drop = FALSE]
    calls = calls[match(samples, calls$sample_id), , drop = FALSE]
    ifelse(calls$assessed, as.integer(calls$value > 0), NA_integer_)
  }
  bcl2 = sv_features$feature[sv_features$gene == "BCL2"][1]
  bcl6 = sv_features$feature[sv_features$gene == "BCL6"][1]

  copy_flat = NULL
  arm = NULL
  has_cn = rep(FALSE, length(samples))
  if(!is.null(seg_data)){
    cn_rows = GAMBLR.utils::get_cnv_features("grch37", feature_set = "LymphGen_model")
    cn_features$direction = cn_rows$direction[match(cn_features$cnv_feature, cn_rows$name)]
    cn_features$arm = cn_rows$region_name[match(cn_features$cnv_feature, cn_rows$name)]
    cn_features$is_arm = cn_rows$type[match(cn_features$cnv_feature, cn_rows$name)] == "arm"
    cn_features = cn_features[!duplicated(cn_features$feature), , drop = FALSE]
    calls = feature_table[feature_table$feature %in% cn_features$feature &
                            feature_table$value > 0, , drop = FALSE]
    has_cn = samples %in% feature_table$sample_id[
      feature_table$feature %in% cn_features$feature & feature_table$assessed]
    calls = calls[calls$assessed, , drop = FALSE]
    f = cn_features[match(calls$feature, cn_features$feature), , drop = FALSE]
    calls$Type = ifelse(f$direction == "loss",
                        ifelse(calls$value >= 2, "HOMDEL", "HETLOSS"),
                        ifelse(calls$value >= 2, "AMP", "GAIN"))
    gene_calls = !f$is_arm
    copy_flat = data.frame(Sample = calls$sample_id[gene_calls],
                           ENTREZ.ID = as.integer(f$entrez[gene_calls]),
                           Type = calls$Type[gene_calls],
                           stringsAsFactors = FALSE)
    arm = data.frame(Sample = calls$sample_id[!gene_calls],
                     Arm = f$arm[!gene_calls],
                     Type = calls$Type[!gene_calls],
                     stringsAsFactors = FALSE)
    copy_flat = copy_flat[order(match(copy_flat$Sample, samples), copy_flat$ENTREZ.ID), , drop = FALSE]
    arm = arm[order(match(arm$Sample, samples), arm$Arm), , drop = FALSE]
    rownames(copy_flat) = NULL
    rownames(arm) = NULL
  }

  sample_annotation = data.frame(Sample.ID = samples,
                                 Copy.Number = as.integer(has_cn),
                                 BCL2.transloc = status(bcl2),
                                 BCL6.transloc = status(bcl6),
                                 stringsAsFactors = FALSE)

  feature_matrix = suppressMessages(GAMBLR.utils::classifier_feature_matrix(
    feature_table, "LymphGen", unassessed_value = NA))

  out = list(mutation_flat = mutation_flat,
             gene_list = data.frame(ENTREZ.ID = gene_ids),
             copy_flat = copy_flat,
             arm = arm,
             sample_annotation = sample_annotation,
             feature_matrix = feature_matrix)
  if(is.null(output_dir)) return(out)

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  file_names = c(mutation_flat = "mutation_flat.tsv", gene_list = "gene_list.txt",
                 copy_flat = "copy_flat.tsv", arm = "arm.tsv",
                 sample_annotation = "sample_annotation.tsv")
  files = character(0)
  for(part in names(file_names)){
    if(is.null(out[[part]])) next
    path = file.path(output_dir, paste0(prefix, "_", file_names[[part]]))
    utils::write.table(out[[part]], path, sep = "\t", quote = FALSE,
                       row.names = FALSE, na = "NA")
    files[part] = path
  }
  attr(out, "files") = files
  invisible(out)
}
