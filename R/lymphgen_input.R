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
#'     (missense, in-frame and nonstop mutations) and Synon (silent and
#'     non-coding mutations of the classes in `synon_classes`, by default
#'     silent and 5' UTR, in a Synon window).
#'     LymphGen counts TRUNC as MUTATION and both as Synon, so each mutation
#'     has one row with its most specific type.
#'   \item Gene list: the Entrez IDs of the genes sequenced. LymphGen refits
#'     its model to these genes, so they must include the unmutated genes.
#'   \item Copy number flat file (Sample, ENTREZ.ID, Type) and arm file
#'     (Sample, Arm, Type): the LymphGen gene and arm copy number features
#'     of the registry, as
#'     HETLOSS or HOMDEL and GAIN or AMP; LymphGen counts HOMDEL as HETLOSS
#'     and AMP as GAIN. Only for samples with copy number. The levels follow
#'     LymphGen's copy number categories: a change of at least half a copy
#'     is HETLOSS or GAIN, and of at least 1.5 copies HOMDEL or AMP, from
#'     the sample's own baseline (from a baseline of 2, HOMDEL is copy
#'     number 0.5 or less and AMP more than 3.5; in a tetraploid tumour a
#'     loss to 2 copies is a HOMDEL), as in LymphGen's training data. Each sample is called on the scale its segments
#'     support ([GAMBLR.utils::cnv_calling_settings]): in absolute copies from
#'     its dominant copy number (`value_scale = "CN"`, `centering = "mode"`)
#'     when they are purity-corrected (a `cn_type` of "absolute", e.g. from
#'     [GAMBLR.utils::calibrate_relative_cn], or else a pipeline such as
#'     battenberg or PureCN in the `pipeline` column), and otherwise as log
#'     ratios centred on the sample's median and converted to copies of a
#'     diploid genome (e.g. gatk_cnv and CNVkit, whose CN is only a rescaled
#'     log ratio). Without those columns, every sample uses CN if
#'     `seg_data` has a CN column, and log ratios otherwise. `cnv_args` overrides the
#'     settings for every sample. Gene features count changes of 30 Mb or
#'     less. For array segments (pipeline "schmitz_array") they are called
#'     as in LymphGen's training data (`gene_rule = "lymphgen"` of
#'     [GAMBLR.utils::segmented_data_to_cnv_status]: one segment of 30 Mb or
#'     less, on its own length, decides each gene); for sequencing (and
#'     segments without a pipeline) a change's length is measured over its
#'     run of neighbouring segments at the same level (`gene_rule =
#'     "fraction"`), since exome segmentation splits broad changes into
#'     pieces under 30 Mb. On the Schmitz exomes linked to LymphGen's
#'     training data, this agrees best with LymphGen's classifications from
#'     the training copy number.
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
#' feature. Without it, the non-coding mutations in `maf_data` are used (e.g.
#' the silent, splice region and targeted region rows of
#' `get_all_coding_ssm(include_silent = TRUE)`).
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
#' @param synon_classes Variant_Classification values of the silent and
#' non-coding mutations that count as Synon (inside the Synon windows); other
#' non-coding mutations are left out of the flat file and the feature matrix.
#' Default `c("Silent", "5'UTR")`, which best matches LymphGen's training
#' data: for the Schmitz exomes linked to it, our in-window Synon calls from
#' these two classes recover 164 of its 199 non-coding Synon features with 99
#' extra, against 170 with 396 extra when intronic, 3' UTR and 5' flank
#' calls also count (most extras are intronic, e.g. in CD83, PIM2, BTG2,
#' PIM1). `c("Silent", "5'UTR", "5'Flank", "Intron", "3'UTR")` reproduces
#' LGenIC's definition.
#' @param synon_regions Where silent and non-coding mutations count as
#' Synon: "registry" (default), the Synon windows of the registry's LymphGen
#' features, which are those LymphGen 2 uses with every gene sequenced; or
#' "gene_tss", the first 4 kb from any active TSS of every gene in the gene
#' list (see [GAMBLR.utils::subset_ssm_by_tss]), for panels, where LymphGen
#' can select other Synon features, and then applies its own windows.
#' @param protein_change_column Column of the MAFs with the protein change.
#' Default NULL uses Protein_Change, or HGVSp_Short if there is none.
#' @param cd79b_location_unknown Set to TRUE to write CD79B's coding
#' mutations (MUTATION and TRUNC rows) with Location -1 ("unknown") in the
#' mutation flat file, so that LymphGen 2 doesn't apply its CD79B position
#' rules to them. Default FALSE. LymphGen 2 drops, in new samples, CD79B
#' mutations past 62007172 and truncations past 62006800 (e.g. the splice
#' site variants beside the ITAM, which are common in MCD). Its training data
#' never had these rules applied: in the first LymphGen release the CD79B
#' threshold was 63007172 (beyond the gene) and the truncation rule compared
#' the wrong column, so neither ever took effect, and LymphGen 2 fixed the
#' rules but kept the same training matrix. In the Schmitz exomes linked to
#' it, the training data have CD79B TRUNC and MUTATION for every such splice
#' site truncation. With TRUE, these mutations count as they did when the
#' model was trained; LymphGen itself is unchanged. Only the mutation flat
#' file is affected, not the feature matrix.
#' @param model The LymphGen model the inputs are for: "original" (default;
#' LymphGen 2 as released) or "gambl_v2" (the GAMBL-derived model, run with
#' `Run_Lymphgen.R --model gambl_v2`, see [run_lymphgen]). "gambl_v2" adds
#' two pseudo-gene features to the mutation flat file, as they were built for
#' its training data:
#' \itemize{
#'   \item NFKBIZ 3' UTR (ENTREZ.ID 1000000006, Type MUTATION): any mutation
#'     in grch37 3:101578215-101578462, whatever its Variant_Classification.
#'     Most coding MAFs leave out 3' UTR mutations, so give them in
#'     `noncoding_maf` (e.g. from `get_ssm_by_region(chromosome = "3",
#'     qstart = 101578215, qend = 101578462)`); without them the feature
#'     counts as unmutated. Exome captures that don't target the UTR still
#'     detect some of these mutations.
#'   \item MYD88 non-L265P (ENTREZ.ID 1000000005, Type MUTATION): MYD88
#'     missense mutations other than L265P.
#' }
#' Both are added to the gene list when their gene (NFKBIZ, MYD88) is in it.
#' "gambl_v3" (`Run_Lymphgen.R --model gambl_v3`) adds one more (the
#' model has no gene-level FOXO1 features; keep FOXO1 in the gene list):
#' \itemize{
#'   \item FOXO1 N-terminal (ENTREZ.ID 1000000007, Type MUTATION): a FOXO1
#'     coding mutation (truncating, missense, in-frame or nonstop class) in
#'     grch37 13:41240271-41240351, the FOXO1_Nterm region of
#'     GAMBLR.utils' curated hotspot regions (the AKT phosphorylation sites,
#'     e.g. R21, S22, T24, and the start codon). Added to the gene list when
#'     FOXO1 is in it.
#' }
#' `Run_Lymphgen.R` drops the model's excluded Synon sub-regions itself. The
#' model has no copy number training data, so `seg_data` is ignored (no
#' copy_flat or arm, Copy.Number 0). The feature matrix is unchanged (the
#' registry's LymphGen features).
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
#' @return A list of data frames, with the model as attribute "lymphgen_model": mutation_flat, gene_list, copy_flat and
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
#'
#' # for the GAMBL-derived model: add the NFKBIZ 3' UTR mutations
#' nfkbiz_utr = get_ssm_by_region(chromosome = "3", qstart = 101578215, qend = 101578462,
#'                                these_samples_metadata = meta, streamlined = FALSE)
#' v2 = lymphgen_input(maf, meta, noncoding_maf = rbind(noncoding, nfkbiz_utr), model = "gambl_v2")
#' run_lymphgen(v2)   # runs Run_Lymphgen.R --model gambl_v2
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
                          synon_classes = c("Silent", "5'UTR"),
                          protein_change_column = NULL,
                          cd79b_location_unknown = FALSE,
                          cnv_args = list(),
                          model = "original",
                          output_dir = NULL,
                          prefix = "lymphgen"){
  model = match.arg(model, c("original", "gambl_v2", "gambl_v3"))
  gambl = model %in% c("gambl_v2", "gambl_v3")
  if(gambl && !is.null(seg_data)){
    message("the ", model, " model has no copy number training data: seg_data is ignored")
    seg_data = NULL
  }
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
  # gambl_v2 pseudo-genes count as sequenced when their gene is
  if(gambl){
    if(4615 %in% gene_ids) gene_ids = c(gene_ids, .gambl_v2$myd88_nonL265P)
    if(64332 %in% gene_ids) gene_ids = c(gene_ids, .gambl_v2$nfkbiz_3utr)   # NFKBIZ is Entrez 64332 (9640, used before, is ZNF592)
  }
  if(model == "gambl_v3" && 2308 %in% gene_ids) gene_ids = c(gene_ids, .gambl_v3$foxo1_nterm)

  # mutations -------------------------------------------------------------
  maf = as.data.frame(maf_data)
  if(!is.null(noncoding_maf)){
    noncoding_maf = as.data.frame(noncoding_maf)
    shared = intersect(colnames(maf), colnames(noncoding_maf))
    maf = rbind(maf[, shared, drop = FALSE], noncoding_maf[, shared, drop = FALSE])
  }
  maf = maf[maf$Tumor_Sample_Barcode %in% samples, , drop = FALSE]
  # no mutations means LymphGen calls every sample Other: say so
  n_mutated = length(unique(maf$Tumor_Sample_Barcode))
  if(n_mutated == 0){
    warning("none of the ", length(samples), " samples has a mutation in maf_data",
            if(!is.null(noncoding_maf)) " or noncoding_maf" else "",
            ", so LymphGen will call them all Other. Check that the MAF's Tumor_Sample_Barcode ",
            "values match these_samples_metadata$sample_id and that its seq_type matches ",
            "(e.g. GAMBLR.open::get_coding_ssm() returns genome samples unless this_seq_type is given)",
            call. = FALSE)
  }else if(n_mutated < length(samples)){
    message(length(samples) - n_mutated, " of ", length(samples), " samples have no mutation in the MAF(s)")
  }
  key_cols =intersect(c("Tumor_Sample_Barcode", "Chromosome", "Start_Position",
                         "Reference_Allele", "Tumor_Seq_Allele2"), colnames(maf))
  maf = maf[!duplicated(maf[, key_cols]), , drop = FALSE]
  maf = GAMBLR.utils::classify_mutations(maf, protein_change_column)
  # gambl_v2 pseudo-gene rows, before non-coding classes are filtered
  extra_flat = if(gambl) .gambl_v2_rows(maf, gene_ids) else NULL
  if(model == "gambl_v3") extra_flat = rbind(extra_flat, .gambl_v3_rows(maf, gene_ids))
  # silent and non-coding mutations count as Synon only in synon_classes; the
  # others are left out of the flat file and the feature matrix
  coding_classes = c("trunc", "missense", "inframe", "nonstop")
  dropped = !maf$mutation_class %in% coding_classes &
    !maf$Variant_Classification %in% synon_classes &
    !(maf$Hugo_Symbol == "MYD88" & maf$Start_Position %in% 38182641)
  maf = maf[!dropped, , drop = FALSE]
  # every sample's mutations, for the feature matrix
  all_maf = maf
  maf$entrez = to_entrez(maf$Hugo_Symbol)
  maf = maf[!is.na(maf$entrez) & maf$entrez %in% gene_ids, , drop = FALSE]
  # coding classes as in the registry's LymphGen MUTATION features; the rest
  # (silent and non-coding mutations of synon_classes) can count as Synon
  coding = maf$mutation_class %in% coding_classes
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
  mutation_flat = rbind(mutation_flat, extra_flat)
  mutation_flat = mutation_flat[order(match(mutation_flat$Sample, samples),
                                      mutation_flat$ENTREZ.ID,
                                      mutation_flat$Location), , drop = FALSE]
  rownames(mutation_flat) = NULL
  # CD79B coding mutations without their position, so LymphGen 2's CD79B
  # position rules (which its training data never had) don't drop them
  if(cd79b_location_unknown){
    cd79b = mutation_flat$ENTREZ.ID %in% 974 & mutation_flat$Type %in% c("MUTATION", "TRUNC")
    mutation_flat$Location[cd79b] = -1L
  }

  # registry features: copy number, translocations and the feature matrix --
  cn_features = registry[registry$data_type == "cnv", , drop = FALSE]
  sv_features = registry[registry$data_type == "sv" &
                           registry$gene %in% c("BCL2", "BCL6"), , drop = FALSE]
  # LymphGen's copy number categories (see details)
  # homdel_max_cn = NULL (relative level 2) as NA, which modifyList keeps
  if("homdel_max_cn" %in% names(cnv_args) && is.null(cnv_args$homdel_max_cn)) cnv_args["homdel_max_cn"] = list(NA_real_)
  user_cnv_args = cnv_args
  has_cn_column = !is.null(seg_data) && "CN" %in% colnames(seg_data)
  cnv_args = utils::modifyList(list(value_scale = if(has_cn_column) "CN" else "log_ratio",
                                    centering = if(has_cn_column) "mode" else "median",
                                    single_thresh = 0.5, high_thresh = 1.5,
                                    gene_rule = "fraction", homdel_max_cn = NA),
                               cnv_args)
  wanted = registry[!duplicated(registry$feature), , drop = FALSE]
  feature_table = suppressMessages(GAMBLR.utils::assemble_classifier_features(
    these_samples_metadata,
    features = wanted[, setdiff(colnames(wanted),
                                c("classifier", "classifier_feature", "group",
                                  "feature_order", "value", "min_level",
                                  "submodels")), drop = FALSE],
    # all_maf has noncoding_maf's mutations too, and without noncoding_maf
    # the non-coding features use its non-coding rows
    maf_data = all_maf,
    seg_data = seg_data,
    sv_data = sv_data,
    sv_from_metadata = sv_from_metadata,
    sv_positive_values = sv_positive_values,
    genome_build = "grch37",
    cnv_args = cnv_args
  ))
  # each sample's copy number features with its own settings, when its
  # segments say what kind of copy number they have (see details)
  groups = .lymphgen_cn_groups(seg_data, samples, user_cnv_args)
  if(!is.null(groups)){
    cn_wanted = wanted[wanted$data_type == "cnv", setdiff(colnames(wanted),
                                                          c("classifier", "classifier_feature", "group",
                                                            "feature_order", "value", "min_level",
                                                            "submodels")), drop = FALSE]
    setting_cols = c("value_scale", "centering", "single_thresh", "high_thresh", "gene_rule", "homdel_max_cn")
    # one key per combination of settings (homdel_max_cn can be NA)
    group_key = do.call(paste, c(groups[, setting_cols], sep = "|"))
    keys = groups[!duplicated(group_key), setting_cols, drop = FALSE]
    for(k in seq_len(nrow(keys))){
      in_group = groups$sample_id[group_key == do.call(paste, c(keys[k, , drop = FALSE], sep = "|"))]
      group_args = utils::modifyList(cnv_args, as.list(keys[k, , drop = FALSE]))
      group_table = suppressMessages(GAMBLR.utils::assemble_classifier_features(
        these_samples_metadata[these_samples_metadata$sample_id %in% in_group, , drop = FALSE],
        features = cn_wanted,
        seg_data = seg_data[seg_data$ID %in% in_group, , drop = FALSE],
        genome_build = "grch37",
        cnv_args = group_args))
      replace = feature_table$feature %in% cn_wanted$feature & feature_table$sample_id %in% in_group
      feature_table = rbind(feature_table[!replace, , drop = FALSE],
                            group_table[, colnames(feature_table), drop = FALSE])
    }
  }
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
  attr(out, "lymphgen_model") = model
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

# The gambl_v2 model's pseudo-genes (as in its training data; see
# README_models.md of the LymphGen repository)
.gambl_v2 = list(nfkbiz_3utr = 1000000006L, nfkbiz_chrom = "3", nfkbiz_start = 101578215, nfkbiz_end = 101578462,
                 myd88_nonL265P = 1000000005L, myd88_l265p_position = 38182641)

# Mutation flat rows of the gambl_v2 pseudo-genes from a classified grch37 MAF:
# NFKBIZ 3' UTR (any mutation in the region) and MYD88 missense other than
# L265P. Only for genes in the gene list.
.gambl_v2_rows = function(maf, gene_ids){
  v = .gambl_v2
  chrom = sub("^chr", "", as.character(maf$Chromosome))
  nfkbiz = chrom == v$nfkbiz_chrom & maf$Start_Position >= v$nfkbiz_start & maf$Start_Position <= v$nfkbiz_end
  myd88 = maf$Hugo_Symbol %in% "MYD88" & maf$Variant_Classification %in% "Missense_Mutation" &
    !(maf$Start_Position %in% v$myd88_l265p_position | maf$protein_change %in% "p.L265P")
  if(!v$nfkbiz_3utr %in% gene_ids) nfkbiz[] = FALSE
  if(!v$myd88_nonL265P %in% gene_ids) myd88[] = FALSE
  message("gambl_v2: NFKBIZ 3' UTR mutations in ", length(unique(maf$Tumor_Sample_Barcode[nfkbiz])),
          " samples; MYD88 non-L265P missense in ", length(unique(maf$Tumor_Sample_Barcode[myd88])), " samples",
          if(!any(nfkbiz)) paste0(". No NFKBIZ 3' UTR mutations: do the MAFs cover grch37 3:", v$nfkbiz_start, "-",
                                  v$nfkbiz_end, "? (see ?lymphgen_input)") else "")
  rows = c(which(nfkbiz), which(myd88))
  data.frame(Sample = maf$Tumor_Sample_Barcode[rows],
             ENTREZ.ID = c(rep(v$nfkbiz_3utr, sum(nfkbiz)), rep(v$myd88_nonL265P, sum(myd88))),
             Type = rep("MUTATION", length(rows)),
             Location = maf$Start_Position[rows],
             stringsAsFactors = FALSE)
}

# The gambl_v3 model's further pseudo-gene: FOXO1 N-terminal, the FOXO1_Nterm
# region of GAMBLR.utils' hotspot_regions.grch37.tsv as used in its training
# data (lymphgen_validation foxo1_nterm_flats.R)
.gambl_v3 = list(foxo1_nterm = 1000000007L, foxo1_chrom = "13", foxo1_start = 41240271, foxo1_end = 41240351)

# Mutation flat rows of FOXO1 N-terminal from a classified grch37 MAF: FOXO1
# coding mutations (trunc, missense, inframe, nonstop classes; i.e. the FOXO1
# MUTATION/TRUNC rows) in the region. Only when the pseudo-gene is in the gene list.
.gambl_v3_rows = function(maf, gene_ids){
  v = .gambl_v3
  chrom = sub("^chr", "", as.character(maf$Chromosome))
  foxo1 = maf$Hugo_Symbol %in% "FOXO1" & chrom == v$foxo1_chrom & maf$Start_Position >= v$foxo1_start &
    maf$Start_Position <= v$foxo1_end & maf$mutation_class %in% c("trunc", "missense", "inframe", "nonstop")
  if(!v$foxo1_nterm %in% gene_ids) foxo1[] = FALSE
  message("gambl_v3: FOXO1 N-terminal mutations in ", length(unique(maf$Tumor_Sample_Barcode[foxo1])), " samples")
  data.frame(Sample = maf$Tumor_Sample_Barcode[foxo1], ENTREZ.ID = rep(v$foxo1_nterm, sum(foxo1)), Type = rep("MUTATION", sum(foxo1)),
             Location = maf$Start_Position[foxo1], stringsAsFactors = FALSE)
}

# Copy number calling settings per sample, from its segments: their cn_type
# column (e.g. from GAMBLR.utils::calibrate_relative_cn) or else their
# pipeline (GAMBLR.utils::cn_pipeline_types), with LymphGen's thresholds;
# a CN-scale sample without CN uses log ratios. NULL when seg_data has
# neither column (every sample then uses the same settings).
.lymphgen_cn_groups = function(seg_data, samples, cnv_args){
  if(is.null(seg_data)) return(NULL)
  segs = as.data.frame(seg_data)
  segs = segs[segs$ID %in% samples, , drop = FALSE]
  if(nrow(segs) == 0 || !any(c("cn_type", "pipeline") %in% colnames(segs))) return(NULL)
  ids = unique(segs$ID)
  pipeline = if("pipeline" %in% colnames(segs)){
    p = tapply(as.character(segs$pipeline), segs$ID, function(x) unique(stats::na.omit(x)))
    several = ids[lengths(p[ids]) > 1]
    if(length(several) > 0){
      stop(length(several), " sample(s) have segments from several pipelines (e.g. ", several[1],
           "); give one segmentation per sample, e.g. from get_cn_segments(source_rule = ...)")
    }
    vapply(p[ids], function(x) if(length(x) == 0) NA_character_ else x, "")
  } else rep(NA_character_, length(ids))
  cn_type = if("cn_type" %in% colnames(segs)){
    t = tapply(as.character(segs$cn_type), segs$ID, function(x) unique(stats::na.omit(x)))
    vapply(t[ids], function(x) if(length(x) == 1) x else NA_character_, "")
  } else rep(NA_character_, length(ids))
  settings = GAMBLR.utils::cnv_calling_settings(ifelse(is.na(pipeline), "unknown", pipeline), "LymphGen",
                                                cnv_args, cn_type = cn_type)
  has_cn = "CN" %in% colnames(segs) &
    ids %in% segs$ID[!is.na(segs$CN)]
  no_cn = settings$value_scale == "CN" & !has_cn
  settings$value_scale[no_cn] = "log_ratio"
  settings$centering[no_cn] = "median"
  data.frame(sample_id = ids, settings[, c("cn_type", "value_scale", "centering", "single_thresh", "high_thresh", "gene_rule",
                                           "homdel_max_cn")],
             stringsAsFactors = FALSE)
}
