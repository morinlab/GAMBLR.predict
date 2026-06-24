#' @title Gene Feature Association Heatmap
#'
#' @description Visualise one or more genes' binary/count feature matrix
#' alongside the coefficients returned by
#' [GAMBLR.utils::select_informative_features]. Each column of the heatmap is
#' one feature (e.g. `MYD88_L265P`, `EZH2_Y641`); each row is a sample with
#' at least one non-zero value across all queried genes. Column annotations
#' below show the glmnet coefficient for every comparison, using a diverging
#' blue–white–red scale centred at zero. When multiple genes are supplied a
#' gene-identity track is prepended to those annotations.
#'
#' @param gene Character scalar or vector of gene symbols. Used as column
#'   prefixes in `all_feats` (e.g. `"MYD88"` matches `MYD88_L265P`,
#'   `MYD88_hotspot`, etc.). When multiple genes are supplied the heatmap
#'   columns are split by gene by default (see `column_split_genes`).
#' @param all_feats Data frame or matrix whose columns are features (`GENE_*`)
#'   and whose rows are samples. Typically the output of
#'   `GAMBLR.utils::summarize_mutation_status()`.
#' @param glm_results Data frame returned by
#'   [GAMBLR.utils::select_informative_features]. Must contain at least the
#'   columns `feature`, `coef`, and `comparison`.
#' @param metadata Optional data frame with sample-level metadata. When
#'   supplied, one or more left-side row annotations are added using the
#'   columns named in `metadata_columns`. Sample matching uses the
#'   `sample_id` column if present, otherwise row names.
#' @param metadata_columns Character vector of column names in `metadata` to
#'   show as row annotations. Colours are assigned automatically via
#'   [GAMBLR.viz::map_metadata_to_colours()]. Default `"pathology"`.
#' @param sort_by_columns Optional character vector of names to sort heatmap
#'   rows (samples) by, in priority order — later entries break ties within
#'   earlier ones. Each name is resolved in order: (1) as a feature matrix
#'   column name — either the full `GENE_suffix` form (e.g. `"EZH2_Y641"`)
#'   or the stripped suffix alone (e.g. `"Y641"`); mutated samples sort first
#'   (descending); (2) as a `metadata` column name (e.g. `"pathology"`,
#'   `"lymphgen"`) — sorted ascending; to control the categorical order
#'   convert the column to an ordered factor with the desired levels before
#'   calling; (3) as a single value present in any metadata column (e.g.
#'   `"FL"`) — samples with that value sort first. Any mix of the three forms
#'   is valid, e.g. `c("pathology", "lymphgen", "EZH2_Y641", "MYD88_L265P")`.
#'   Overrides `cluster_rows`. Default `NULL` (use clustering).
#' @param strip_gene_prefix Logical. If `TRUE` (default), removes the leading
#'   `GENE_` prefix from column labels. Gene identity is still conveyed via
#'   the per-block titles when `column_split_genes = TRUE`.
#' @param drop_uninformative Logical. If `TRUE` (default), features with a
#'   zero coefficient in every comparison are dropped and a warning is issued.
#'   Set to `FALSE` to retain all features.
#' @param sort_features_by Optional character vector of comparison names
#'   (after `"vs rest"` trimming, e.g. `c("DLBCL", "FL")`). Columns are
#'   sorted by descending coefficient for the first group, ties broken by the
#'   second, and so on. When `column_split_genes = TRUE` the sort is applied
#'   independently within each gene block; otherwise it is global. Overrides
#'   `cluster_columns`. Default `NULL` (use clustering).
#' @param column_split_genes Logical. When multiple genes are supplied,
#'   whether to use ComplexHeatmap's `column_split` to visually group features
#'   by gene. Defaults to `TRUE` for multiple genes, `FALSE` for one.
#' @param show_annotation_legend Logical. If `TRUE` (default), a single colour
#'   legend for the coefficient scale is shown (attached to the first
#'   comparison track). Set to `FALSE` to suppress it entirely.
#' @param show_heatmap_legend Logical. Whether to show the body colour scale
#'   legend. Default `FALSE` — suppressed because values are typically binary
#'   and the scale adds little information.
#' @param heatmap_title Optional character string used as the heatmap column
#'   title. Defaults to the gene name(s) joined by `"/"`, which provides a
#'   label when the body legend is suppressed and `strip_gene_prefix = TRUE`
#'   would otherwise leave the genes unidentified.
#' @param annotation_bar_height Height of each annotation track in millimetres,
#'   applied to both the coefficient tracks and the metadata tracks via
#'   `simple_anno_size`. Default `3`.
#' @param annotation_name_fontsize Font size (pt) for annotation track name
#'   labels. Default `8`.
#' @param rotate Logical. If `TRUE`, transposes the heatmap so features are
#'   rows and samples are columns. The coefficient annotation moves to the
#'   right side and the metadata annotation moves to the bottom. Default
#'   `FALSE`.
#' @param hide_sample_id Logical. If `TRUE` (default), sample identifiers are
#'   not printed on the heatmap. When `rotate = FALSE` this suppresses row
#'   names; when `rotate = TRUE` it suppresses column names.
#' @param show_feature_names Logical. Whether to display feature names on the
#'   heatmap. When `rotate = FALSE` features are columns (`show_column_names`);
#'   when `rotate = TRUE` features are rows (`show_row_names`). Default `TRUE`.
#' @param feature_order Optional character vector of feature names (using the
#'   full `GENE_suffix` form, before any prefix stripping) specifying the exact
#'   column order. Features absent from the vector are dropped. Overrides both
#'   `sort_features_by` and `cluster_columns`. Default `NULL`.
#' @param gene_order Optional character vector of gene names controlling the
#'   left-to-right order of gene blocks when `column_split_genes = TRUE`.
#'   Must be a permutation of the supplied `gene` vector. Defaults to the
#'   order genes are supplied in `gene`.
#' @param cluster_columns Logical. Cluster heatmap columns. Default `TRUE`.
#' @param cluster_rows Logical. Cluster heatmap rows. Default `TRUE`.
#' @param order_by_metadata Logical. When `TRUE` and `metadata` is supplied,
#'   rows (samples) are ordered to match the row order of `metadata` rather
#'   than being clustered or sorted via `sort_by_columns`. Row clustering is
#'   disabled automatically; `sort_by_columns` is ignored. Feature ordering
#'   (`sort_features_by`, `feature_order`, `cluster_columns`) is unaffected.
#'   Default `FALSE`.
#' @param ... Additional arguments passed directly to
#'   [ComplexHeatmap::Heatmap()], e.g. `rect_gp`.
#'
#' @return A [ComplexHeatmap::Heatmap-class] object.
#'
#' @import ComplexHeatmap
#' @import dplyr
#' @export
#'
#' @examples
#' \dontrun{
#' glm_res <- select_informative_features(feature_df, meta)
#'
#' # single gene
#' gene_feat_heatmap("MYD88", feature_df, glm_res)
#'
#' # multiple genes, features split by gene
#' gene_feat_heatmap(c("MYD88", "EZH2"), feature_df, glm_res)
#' }
gene_feat_heatmap <- function(
  gene,
  all_feats,
  glm_results,
  metadata               = NULL,
  metadata_columns       = "pathology",
  sort_by_columns        = NULL,
  strip_gene_prefix      = NULL,
  drop_uninformative     = TRUE,
  sort_features_by       = NULL,
  column_split_genes     = NULL,
  show_annotation_legend = TRUE,
  show_heatmap_legend    = FALSE,
  heatmap_title          = NULL,
  annotation_bar_height    = 3,
  annotation_name_fontsize = 8,
  rotate                   = FALSE,
  hide_sample_id           = TRUE,
  show_feature_names       = TRUE,
  feature_order          = NULL,
  gene_order             = NULL,
  cluster_columns        = TRUE,
  cluster_rows           = TRUE,
  order_by_metadata      = FALSE,
  ...
) {
  genes      <- unique(gene)
  multi_gene <- length(genes) > 1L

  # resolve defaults that depend on single vs multi-gene mode
  if (is.null(strip_gene_prefix))  strip_gene_prefix  <- TRUE
  if (is.null(column_split_genes)) column_split_genes <- multi_gene

  # ── build feature matrix, one gene at a time ──────────────────────────────
  gene_col_list <- lapply(genes, function(g) {
    cols <- dplyr::select(all_feats, dplyr::starts_with(paste0(g, "_")))
    if (ncol(cols) == 0L) {
      warning("No features found for gene '", g, "'; skipping.")
      return(NULL)
    }
    cols
  })
  names(gene_col_list) <- genes
  gene_col_list        <- Filter(Negate(is.null), gene_col_list)
  if (length(gene_col_list) == 0L) {
    stop("No features found for any of the supplied gene(s).")
  }

  # col_gene[i] = which gene column i belongs to; kept in sync throughout
  col_gene  <- rep(
    names(gene_col_list),
    times = vapply(gene_col_list, ncol, integer(1L))
  )
  gene_cols <- dplyr::bind_cols(gene_col_list)

  sums     <- rowSums(gene_cols, na.rm = TRUE)
  mut_rows <- as.matrix(gene_cols[sums > 0, , drop = FALSE])
  if (nrow(mut_rows) == 0L) {
    stop("No mutated samples found for ",
         paste(genes, collapse = "/"), ".")
  }

  # ── heatmap body colour (white → navy) ───────────────────────────────────
  body_max <- max(mut_rows, na.rm = TRUE)
  if (body_max == 0) body_max <- 1L
  body_col <- circlize::colorRamp2(c(0, body_max), c("white", "#1F4E79"))

  # ── coefficient matrix (all genes × all comparisons) ─────────────────────
  gene_results <- dplyr::filter(glm_results,
                                .data$feature %in% colnames(gene_cols))

  coef_anno <- NULL
  if (nrow(gene_results) > 0 &&
        all(c("feature", "coef", "comparison") %in% colnames(gene_results))) {

    coef_wide <- tidyr::pivot_wider(
      gene_results,
      id_cols     = "feature",
      names_from  = "comparison",
      values_from = "coef",
      values_fill = 0
    )
    coef_mat           <- as.matrix(dplyr::select(coef_wide, -"feature"))
    rownames(coef_mat) <- coef_wide$feature

    # pad features absent from glm output with coef = 0
    missing_feats <- setdiff(colnames(mut_rows), rownames(coef_mat))
    if (length(missing_feats) > 0L) {
      filler <- matrix(
        0,
        nrow     = length(missing_feats),
        ncol     = ncol(coef_mat),
        dimnames = list(missing_feats, colnames(coef_mat))
      )
      coef_mat <- rbind(coef_mat, filler)
    }
    coef_mat <- coef_mat[colnames(mut_rows), , drop = FALSE]

    # drop uninformative features; keep col_gene in sync via logical mask
    if (drop_uninformative) {
      keep_mask <- rowSums(abs(coef_mat)) > 0
      zero_feats <- rownames(coef_mat)[!keep_mask]
      if (length(zero_feats) > 0L) {
        warning(
          length(zero_feats),
          " feature(s) with zero coefficient in all comparisons dropped: ",
          paste(zero_feats, collapse = ", ")
        )
        mut_rows <- mut_rows[, keep_mask, drop = FALSE]
        coef_mat <- coef_mat[keep_mask, , drop = FALSE]
        col_gene <- col_gene[keep_mask]
      }
    }

    # trim "vs rest" from comparison labels
    colnames(coef_mat) <- sub(" vs rest$", "", colnames(coef_mat))

    # column sort — within gene blocks when splitting, otherwise global
    if (!is.null(sort_features_by)) {
      missing_groups <- setdiff(sort_features_by, colnames(coef_mat))
      if (length(missing_groups) > 0L) {
        warning(
          "Groups in 'sort_features_by' not found in comparisons ",
          "and ignored: ",
          paste(missing_groups, collapse = ", ")
        )
      }
      valid_groups <- intersect(sort_features_by, colnames(coef_mat))
      if (length(valid_groups) > 0L) {
        if (column_split_genes && multi_gene) {
          # sort within each gene block, preserving block order
          blocks    <- split(seq_len(ncol(mut_rows)), col_gene)
          col_order <- unlist(lapply(genes, function(g) {
            idx  <- blocks[[g]]
            if (is.null(idx) || length(idx) == 0L) return(integer(0L))
            keys <- lapply(valid_groups, function(grp) -coef_mat[idx, grp])
            idx[do.call(order, keys)]
          }), use.names = FALSE)
        } else {
          keys      <- lapply(valid_groups, function(g) -coef_mat[, g])
          col_order <- do.call(order, keys)
        }
        mut_rows        <- mut_rows[, col_order, drop = FALSE]
        coef_mat        <- coef_mat[col_order, , drop = FALSE]
        col_gene        <- col_gene[col_order]
        cluster_columns <- FALSE
      }
    }

    # explicit feature order — overrides sort_features_by and clustering
    if (!is.null(feature_order)) {
      unknown <- setdiff(feature_order, colnames(mut_rows))
      if (length(unknown) > 0L) {
        warning("Features in 'feature_order' not found and ignored: ",
                paste(unknown, collapse = ", "))
      }
      keep <- intersect(feature_order, colnames(mut_rows))
      if (length(keep) > 0L) {
        idx             <- match(keep, colnames(mut_rows))
        mut_rows        <- mut_rows[, idx, drop = FALSE]
        coef_mat        <- coef_mat[idx, , drop = FALSE]
        col_gene        <- col_gene[idx]
        cluster_columns <- FALSE
      }
    }

    # shared diverging colour scale across all comparisons
    coef_lim <- max(abs(coef_mat), na.rm = TRUE)
    if (coef_lim == 0) coef_lim <- 1
    coef_col <- circlize::colorRamp2(
      c(-coef_lim, 0, coef_lim),
      c("#2166AC", "white", "#D6604D")
    )

    # build coefficient annotation tracks
    coef_data     <- lapply(seq_len(ncol(coef_mat)),
                            function(j) coef_mat[, j])
    names(coef_data) <- colnames(coef_mat)

    coef_col_list <- stats::setNames(
      rep(list(coef_col), ncol(coef_mat)),
      colnames(coef_mat)
    )
    show_leg        <- rep(FALSE, ncol(coef_mat))
    if (show_annotation_legend) show_leg[1] <- TRUE
    names(show_leg) <- colnames(coef_mat)

    if (multi_gene) {
      # gene identity track: one colour per gene
      gambl_cols <- GAMBLR.helpers::get_gambl_colours()
      have_cols  <- intersect(genes, names(gambl_cols))
      need_cols  <- setdiff(genes, names(gambl_cols))
      gene_pal   <- gambl_cols[have_cols]
      if (length(need_cols) > 0L) {
        fallback <- grDevices::hcl.colors(length(need_cols), "Dark 2")
        gene_pal <- c(gene_pal,
                      stats::setNames(fallback, need_cols))
      }
      gene_pal <- gene_pal[genes]  # consistent gene order

      gene_track    <- factor(col_gene, levels = genes)
      anno_data     <- c(list(Gene = gene_track), coef_data)
      col_arg       <- c(list(Gene = gene_pal), coef_col_list)
      show_leg_full <- c(Gene = TRUE, show_leg)
    } else {
      anno_data     <- coef_data
      col_arg       <- coef_col_list
      show_leg_full <- show_leg
    }

    # override the legend title on the first visible coef track to "Coef"
    coef_legend_param <- stats::setNames(
      list(list(title = "Coef")),
      colnames(coef_mat)[1]
    )

    coef_anno_fn <- if (rotate) {
      ComplexHeatmap::rowAnnotation
    } else {
      ComplexHeatmap::columnAnnotation
    }
    coef_anno <- do.call(
      coef_anno_fn,
      c(anno_data, list(
        col                     = col_arg,
        show_legend             = show_leg_full,
        annotation_legend_param = coef_legend_param,
        simple_anno_size        = grid::unit(annotation_bar_height, "mm"),
        annotation_name_gp   = grid::gpar(fontsize = annotation_name_fontsize),
        annotation_name_side    = if (rotate) "top" else "left"
      ))
    )

  } else {
    warning(
      "No glm_results entries found for gene(s) '",
      paste(genes, collapse = "', '"),
      "'; drawing heatmap without coefficient annotation."
    )
    if (!is.null(feature_order)) {
      unknown <- setdiff(feature_order, colnames(mut_rows))
      if (length(unknown) > 0L) {
        warning("Features in 'feature_order' not found and ignored: ",
                paste(unknown, collapse = ", "))
      }
      keep <- intersect(feature_order, colnames(mut_rows))
      if (length(keep) > 0L) {
        idx             <- match(keep, colnames(mut_rows))
        mut_rows        <- mut_rows[, idx, drop = FALSE]
        col_gene        <- col_gene[idx]
        cluster_columns <- FALSE
      }
    }
  }

  # ── resolve metadata sample index ─────────────────────────────────────────
  meta_idx <- NULL
  if (!is.null(metadata)) {
    id_col   <- if ("sample_id" %in% colnames(metadata)) "sample_id" else NULL
    meta_idx <- if (!is.null(id_col)) {
      match(rownames(mut_rows), metadata[[id_col]])
    } else {
      match(rownames(mut_rows), rownames(metadata))
    }

    # drop samples absent from metadata
    in_meta <- !is.na(meta_idx)
    if (any(!in_meta)) {
      message(sum(!in_meta),
              " sample(s) not found in metadata and dropped from plot.")
      mut_rows <- mut_rows[in_meta, , drop = FALSE]
      meta_idx <- meta_idx[in_meta]
    }
  }

  # ── metadata-order mode: reorder rows to follow metadata row sequence ────
  if (order_by_metadata && !is.null(metadata)) {
    id_col   <- if ("sample_id" %in% colnames(metadata)) "sample_id" else NULL
    meta_ids <- if (!is.null(id_col)) metadata[[id_col]] else rownames(metadata)
    row_ids  <- rownames(mut_rows)
    # walk metadata in order, keeping only samples present in mut_rows
    ordered_ids <- meta_ids[meta_ids %in% row_ids]
    row_order   <- match(ordered_ids, row_ids)
    row_order   <- row_order[!is.na(row_order)]
    mut_rows    <- mut_rows[row_order, , drop = FALSE]
    if (!is.null(meta_idx)) meta_idx <- meta_idx[row_order]
    cluster_rows <- FALSE
  }

  # ── optional row sort: feature columns and/or metadata ────────────────────
  if (!order_by_metadata && !is.null(sort_by_columns)) {
    sort_keys  <- list()
    unresolved <- character(0)
    # precompute stripped column names so users can pass e.g. "Y641" for "EZH2_Y641"
    stripped_cols <- colnames(mut_rows)
    for (g in genes) stripped_cols <- sub(paste0("^", g, "_"), "", stripped_cols)
    for (col in sort_by_columns) {
      if (col %in% colnames(mut_rows)) {
        sort_keys[[col]] <- -mut_rows[, col]
      } else if (col %in% stripped_cols) {
        hits <- which(stripped_cols == col)
        if (length(hits) > 1L) {
          warning("'", col, "' matches ", length(hits),
                  " features after prefix stripping; using '",
                  colnames(mut_rows)[hits[1L]], "'.")
        }
        sort_keys[[col]] <- -mut_rows[, hits[1L]]
      } else if (!is.null(meta_idx) && col %in% colnames(metadata)) {
        sort_keys[[col]] <- metadata[[col]][meta_idx]
      } else if (!is.null(meta_idx)) {
        found <- FALSE
        for (mcol in colnames(metadata)) {
          vals <- as.character(metadata[[mcol]][meta_idx])
          if (col %in% vals) {
            sort_keys[[col]] <- -(vals == col)
            found <- TRUE
            break
          }
        }
        if (!found) unresolved <- c(unresolved, col)
      } else {
        unresolved <- c(unresolved, col)
      }
    }
    if (length(unresolved) > 0L) {
      warning(
        "Names in 'sort_by_columns' not found in feature matrix, ",
        "metadata column names, or metadata values and ignored: ",
        paste(unresolved, collapse = ", ")
      )
    }
    if (length(sort_keys) > 0L) {
      row_order    <- do.call(order, sort_keys)
      mut_rows     <- mut_rows[row_order, , drop = FALSE]
      if (!is.null(meta_idx)) meta_idx <- meta_idx[row_order]
      cluster_rows <- FALSE
    }
  }

  # ── left-side metadata row annotation ─────────────────────────────────────
  meta_anno <- NULL
  if (!is.null(metadata)) {
    missing_cols <- setdiff(metadata_columns, colnames(metadata))
    if (length(missing_cols) > 0L) {
      warning("Columns not found in metadata and ignored: ",
              paste(missing_cols, collapse = ", "))
      metadata_columns <- intersect(metadata_columns, colnames(metadata))
    }
    if (length(metadata_columns) > 0L) {
      meta_aligned <- as.data.frame(
        metadata[meta_idx, metadata_columns, drop = FALSE]
      )
      rownames(meta_aligned) <- rownames(mut_rows)
      colour_list <- GAMBLR.viz::map_metadata_to_colours(
        metadataColumns        = metadata_columns,
        these_samples_metadata = meta_aligned
      )
      meta_anno_fn <- if (rotate) {
        ComplexHeatmap::columnAnnotation
      } else {
        ComplexHeatmap::rowAnnotation
      }
      meta_anno <- meta_anno_fn(
        df                   = meta_aligned,
        col                  = colour_list,
        simple_anno_size     = grid::unit(annotation_bar_height, "mm"),
        annotation_name_gp   = grid::gpar(fontsize = annotation_name_fontsize),
        annotation_name_side = if (rotate) "left" else "top"
      )
    }
  }

  # ── column split factor (multi-gene only) ─────────────────────────────────
  split_levels <- if (!is.null(gene_order)) gene_order else genes
  col_split <- if (column_split_genes && multi_gene) {
    factor(col_gene, levels = split_levels)
  } else {
    NULL
  }

  # ── strip per-gene prefix from column labels (after all sorting) ──────────
  if (strip_gene_prefix) {
    new_names <- colnames(mut_rows)
    for (g in genes) {
      new_names <- sub(paste0("^", g, "_"), "", new_names)
    }
    colnames(mut_rows) <- new_names
  }

  # ── assemble heatmap ──────────────────────────────────────────────────────
  if (is.null(heatmap_title)) heatmap_title <- paste(genes, collapse = "/")

  # when gene splitting is active, title each block with its gene name;
  # otherwise use heatmap_title for the single overall title
  split_titles <- if (column_split_genes && multi_gene) split_levels else heatmap_title

  if (rotate) {
    ComplexHeatmap::Heatmap(
      t(mut_rows),
      name                = paste(genes, collapse = "/"),
      col                 = body_col,
      row_title           = split_titles,
      show_heatmap_legend = show_heatmap_legend,
      show_row_names      = show_feature_names,
      show_column_names   = !hide_sample_id,
      cluster_rows        = cluster_columns,
      cluster_columns     = cluster_rows,
      row_split           = col_split,
      right_annotation    = coef_anno,
      bottom_annotation   = meta_anno,
      ...
    )
  } else {
    ComplexHeatmap::Heatmap(
      mut_rows,
      name                = paste(genes, collapse = "/"),
      col                 = body_col,
      column_title        = split_titles,
      show_heatmap_legend = show_heatmap_legend,
      show_row_names      = !hide_sample_id,
      show_column_names   = show_feature_names,
      cluster_columns     = cluster_columns,
      cluster_rows        = cluster_rows,
      column_split        = col_split,
      left_annotation     = meta_anno,
      bottom_annotation   = coef_anno,
      ...
    )
  }
}
