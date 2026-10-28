#' @title Build a DLBClass genetic sample matrix (GSM)
#'
#' @description Build the input for the DLBClass classifier (the "GSM": 163
#' mutation, structural variant and copy number features) from GAMBLR
#' pipeline outputs, using the existing GAMBLR feature functions:
#' [GAMBLR.utils::assemble_genetic_features] for mutations and SVs and
#' [GAMBLR.utils::segmented_data_to_cnv_status] for copy number. Feature
#' names, encodings and feature order follow DLBClass (the DLBclass-tool
#' scripts maf2gsm.py, sv2gsm.py, seg2gsm.py and combine2gsm.py).
#'
#' @details Encodings, as in DLBClass:
#' \itemize{
#'   \item Mutations (one column per gene): as in maf2gsm.py, only
#'     mutations with a protein change (`protein_change_column`) or a
#'     splice site / start codon class are used. Those other than Silent are
#'     non-silent (2, including splice region variants with a protein
#'     change); otherwise a Silent mutation gives 1; otherwise 0. MYD88 is also
#'     split into `MYD88.L265P` (2 for p.L265P) and `MYD88.OTHER` (2 for any
#'     other non-silent mutation, 1 for a silent one). `maf_data` should
#'     include silent mutations (e.g. `get_coding_ssm(include_silent =
#'     TRUE)`).
#'   \item Structural variants (`SV.BCL2`, `SV.BCL6`, `SV.MYC`): 3 when the
#'     metadata column named in `sv_from_metadata` is positive, otherwise 0
#'     (including unknown).
#'   \item Copy number (`X17P.DEL` etc.): the DLBClass arm-level and focal
#'     features called by [GAMBLR.utils::segmented_data_to_cnv_status] (0, 1
#'     or 2). Without `seg_data` they are all 0, with a warning.
#' }
#' Features with no events in the cohort are included as all-zero columns,
#' and every sample in `these_samples_metadata` is included.
#'
#' @param maf_data MAF data frame of coding (and silent) mutations for the
#' samples, e.g. from `get_coding_ssm(..., include_silent = TRUE)`.
#' @param these_samples_metadata Metadata with one row per sample_id; defines
#' the samples (mRNA-seq samples are dropped, with a warning) and supplies the
#' SV status columns.
#' @param seg_data Optional segmented copy number data (columns ID, chrom,
#' start, end, log.ratio) for the copy number features.
#' @param genome_build Genome build of `seg_data`. Default "grch37" (the only
#' build [GAMBLR.utils::segmented_data_to_cnv_status] currently supports).
#' @param protein_change_column Column of `maf_data` with the protein change
#' (e.g. "p.L265P"). Default NULL uses Protein_Change, or HGVSp_Short (as in
#' GAMBL MAFs) if there is no Protein_Change column.
#' @param sv_from_metadata Named vector mapping the SV features to metadata
#' columns. Default `c(BCL2 = "bcl2_ba", BCL6 = "bcl6_ba", MYC = "myc_ba")`.
#' @param sv_positive_values Metadata values that indicate an SV. Default
#' "POS".
#' @param cnv_args Optional named list of further arguments for
#' [GAMBLR.utils::segmented_data_to_cnv_status], e.g. `blacklist`.
#' @param output_format "matrix" (default; samples as rows, features as
#' columns) or "gsm" (DLBClass layout: a classifier_name column, then one
#' column per sample, as written by DLBclass-tool).
#'
#' @return A data frame in the requested format, with the 163 DLBClass
#' features in DLBClass's order.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' meta = get_gambl_metadata() %>% dplyr::filter(pathology == "DLBCL")
#' maf = get_coding_ssm(these_samples_metadata = meta, include_silent = TRUE)
#' segs = get_cn_segments(these = meta, projection = "grch37")
#' gsm = dlbclass_gsm(maf, meta, seg_data = segs)
#'
#' # write the DLBClass layout for DLBclass-tool
#' readr::write_tsv(dlbclass_gsm(maf, meta, seg_data = segs,
#'                               output_format = "gsm"),
#'                  "cohort.GSM.tsv")
#' }
dlbclass_gsm = function(maf_data,
                        these_samples_metadata,
                        seg_data = NULL,
                        genome_build = "grch37",
                        protein_change_column = NULL,
                        sv_from_metadata = c(BCL2 = "bcl2_ba",
                                             BCL6 = "bcl6_ba",
                                             MYC = "myc_ba"),
                        sv_positive_values = "POS",
                        cnv_args = list(),
                        output_format = "matrix"){
  if(!output_format %in% c("matrix", "gsm")){
    stop("output_format must be \"matrix\" or \"gsm\"")
  }
  if(!"sample_id" %in% colnames(these_samples_metadata)){
    stop("these_samples_metadata needs a sample_id column")
  }
  feature_order = readLines(system.file("extdata", "dlbclass_feature_order.txt",
                                        package = "GAMBLR.predict"))
  these_samples_metadata = GAMBLR.utils::drop_mrna_samples(these_samples_metadata,
                                                           "dlbclass_gsm")
  samples = unique(as.character(these_samples_metadata$sample_id))
  # DLBClass gene names, in GAMBLR's spelling (HLA-A rather than HLA.A)
  ssm_features = feature_order[!grepl("^X|^SV[.]", feature_order)]
  genes = setdiff(ssm_features, c("MYD88.L265P", "MYD88.OTHER"))
  genes = sub("^HLA[.]", "HLA-", genes)

  # mutations, selected as in maf2gsm.py: a protein change or a splice site /
  # start codon class
  maf_data = as.data.frame(maf_data)
  if(is.null(protein_change_column)){
    protein_change_column = intersect(c("Protein_Change", "HGVSp_Short"),
                                      colnames(maf_data))[1]
  }
  if(is.na(protein_change_column) ||
     !protein_change_column %in% colnames(maf_data)){
    stop("maf_data needs a protein change column (Protein_Change or HGVSp_Short), or set protein_change_column")
  }
  protein_change = as.character(maf_data[[protein_change_column]])
  functional_classes = c("Splice_Site", "De_novo_Start_InFrame",
                         "De_novo_Start_OutOfFrame", "Start_Codon_Del")
  used = (!is.na(protein_change) & protein_change != "") |
    maf_data$Variant_Classification %in% functional_classes
  maf_data = maf_data[used & maf_data$Tumor_Sample_Barcode %in% samples &
                        maf_data$Hugo_Symbol %in% genes, , drop = FALSE]
  if(!any(maf_data$Variant_Classification == "Silent")){
    message("maf_data has no silent mutations; use include_silent = TRUE to ",
            "encode them as DLBClass does (1)")
  }
  # DLBClass counts every selected mutation other than Silent as non-silent,
  # while assemble_genetic_features counts only vc_nonSynonymous classes as
  # coding, so relabel the others (e.g. splice region variants with a protein
  # change) to a non-synonymous class
  vc = as.character(maf_data$Variant_Classification)
  relabel = vc != "Silent" & !vc %in% GAMBLR.helpers::vc_nonSynonymous
  vc[relabel & vc == "Splice_Region"] = "Splice_Site"
  vc[relabel & vc %in% c("De_novo_Start_InFrame", "De_novo_Start_OutOfFrame",
                         "Start_Codon_Del")] = "Translation_Start_Site"
  vc[relabel & !vc %in% GAMBLR.helpers::vc_nonSynonymous] = "Missense_Mutation"
  maf_data$Variant_Classification = vc

  # non-silent mutations are coding (2) and silent ones non-coding (1), and
  # each gene's columns are combined by their maximum, as DLBClass encodes
  # them. MYD88 is split into its hotspot (L265P) and other mutations.
  status = NULL
  if(nrow(maf_data) > 0){
    # without MYD88 mutations, annotate_curated_drivers warns that it has
    # nothing to annotate, which is expected here
    status = withCallingHandlers(GAMBLR.utils::assemble_genetic_features(
      maf_data,
      genes_coding = genes,
      genes_noncoding = genes,
      genes_hotspot = "MYD88",
      genes_driver = character(0),
      encoding_policy = list(coding = 2, driver = 2, noncoding = 1, sv = 0),
      these_samples_metadata = these_samples_metadata,
      sv_from_metadata = NULL
    ), warning = function(w){
      if(grepl("No mutations found in genes of interest", conditionMessage(w))){
        invokeRestart("muffleWarning")
      }
    })
  }

  gsm = matrix(0, nrow = length(samples), ncol = length(feature_order),
               dimnames = list(samples, feature_order))
  if(!is.null(status)){
    status = as.matrix(status)
    rows = intersect(rownames(status), samples)
    for(col in colnames(status)){
      # DLBClass names: HLA.A rather than HLA-A
      target = sub("^HLA-", "HLA.", col)
      if(target %in% feature_order){
        gsm[rows, target] = status[rows, col]
      }
    }
    # MYD88: any mutation (MYD88), p.L265P (MYD88.L265P, the hotspot) and any
    # other mutation, including silent ones (MYD88.OTHER)
    column_max = function(cols){
      if(length(cols) == 0) return(rep(0, length(rows)))
      apply(status[rows, cols, drop = FALSE], 1, max)
    }
    myd88_any = grep("^MYD88(_|$)", colnames(status), value = TRUE)
    gsm[rows, "MYD88"] = column_max(myd88_any)
    gsm[rows, "MYD88.OTHER"] = column_max(setdiff(myd88_any, "MYD88_hotspot"))
    gsm[rows, "MYD88.L265P"] = column_max(intersect("MYD88_hotspot", myd88_any))
  }

  # structural variants from metadata (3 when positive)
  if(!is.null(sv_from_metadata) && length(sv_from_metadata) > 0){
    meta_rows = match(samples, as.character(these_samples_metadata$sample_id))
    for(gene in names(sv_from_metadata)){
      column = sv_from_metadata[[gene]]
      if(!column %in% colnames(these_samples_metadata)){
        warning("SV metadata column ", column, " not found; SV.", gene,
                " is 0 for all samples")
        next
      }
      positive = as.character(these_samples_metadata[[column]][meta_rows]) %in%
        sv_positive_values
      gsm[, paste0("SV.", gene)] = ifelse(positive, 3, 0)
    }
  }

  # copy number features
  if(!is.null(seg_data)){
    cnv_call = utils::modifyList(cnv_args, list(
      seg_data = seg_data,
      cnv_features = GAMBLR.utils::get_cnv_features(genome_build,
                                                    feature_set = "DLBClass"),
      these_samples_metadata = these_samples_metadata,
      genome_build = genome_build,
      output_format = "gsm"
    ))
    cnv = do.call(GAMBLR.utils::segmented_data_to_cnv_status, cnv_call)
    cnv_names = paste0("X", cnv$classifier_name)
    unknown = setdiff(cnv_names, feature_order)
    if(length(unknown) > 0){
      warning("copy number features not in the DLBClass feature order: ",
              paste(unknown, collapse = ", "))
    }
    cnv_values = as.matrix(cnv[, setdiff(colnames(cnv), "classifier_name"),
                               drop = FALSE])
    rownames(cnv_values) = cnv_names
    cnv_samples = intersect(colnames(cnv_values), samples)
    keep = intersect(cnv_names, feature_order)
    gsm[cnv_samples, keep] = t(cnv_values[keep, cnv_samples, drop = FALSE])
  }else{
    warning("no seg_data: the DLBClass copy number features are all 0")
  }

  if(output_format == "gsm"){
    out = data.frame(classifier_name = feature_order, t(gsm),
                     check.names = FALSE)
    rownames(out) = NULL
    return(out)
  }
  as.data.frame(gsm)
}
