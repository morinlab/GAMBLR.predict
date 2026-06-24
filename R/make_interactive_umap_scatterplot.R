#' @title Interactive UMAP Scatterplot With Feature Hover
#'
#' @description Generates an interactive plotly scatterplot from the output of
#' `make_and_annotate_umap()`, mirroring the aesthetics of
#' [GAMBLR.predict::make_umap_scatterplot]. Hovering over any point displays
#' all non-zero (mutated) features for that sample, drawn from the feature
#' matrix stored in the result object.
#'
#' @param umap_result The list returned by `make_and_annotate_umap()`. Must
#'   contain a `$df` element (UMAP coordinates and metadata, with columns `V1`
#'   and `V2`) and a `$features` element (the sample-by-feature matrix used to
#'   build the UMAP, with sample IDs as row names).
#' @param sample_id_column Name of the column in `umap_result$df` that holds
#'   sample identifiers, used to match rows in `$features`. Default
#'   `"sample_id"`.
#' @param drop_composite Logical. If `TRUE`, removes samples whose
#'   `colour_by` value contains `"COMP"`. Default `TRUE`.
#' @param colour_by Name of the column in `umap_result$df` used to colour
#'   points. Default `"lymphgen"`.
#' @param drop_other Logical. If `TRUE`, removes samples with `colour_by`
#'   equal to `"Other"` or `"NOS"`. Default `FALSE`.
#' @param high_confidence Logical. If `TRUE`, retains only samples with
#'   `Confidence > 0.7`. Default `FALSE`.
#' @param custom_colours Optional named character vector mapping group labels
#'   to hex colours. Defaults to [GAMBLR.helpers::get_gambl_colours()].
#' @param title Optional plot title string.
#' @param alpha Point opacity (0–1). Default `0.8`.
#' @param point_size Point diameter in pixels. Default `5`.
#' @param max_features Maximum number of feature names to show per sample in
#'   the hover tooltip. Features beyond this limit are summarised as
#'   `"... +N more"`. Default `30`.
#'
#' @return A `plotly` htmlwidget.
#'
#' @import dplyr
#' @export
#'
#' @examples
#' \dontrun{
#' my_umap <- make_and_annotate_umap(feature_df, metadata)
#'
#' make_interactive_umap_scatterplot(
#'   umap_result = my_umap,
#'   colour_by   = "lymphgen",
#'   drop_other  = TRUE
#' )
#' }
make_interactive_umap_scatterplot <- function(
  umap_result,
  sample_id_column = "sample_id",
  drop_composite   = TRUE,
  colour_by        = "lymphgen",
  drop_other       = FALSE,
  high_confidence  = FALSE,
  custom_colours,
  title            = NULL,
  alpha            = 0.8,
  point_size       = 5,
  max_features     = 30
) {
  if (!is.list(umap_result) || !all(c("df", "features") %in% names(umap_result))) {
    stop("'umap_result' must be the list returned by make_and_annotate_umap(), ",
         "containing '$df' and '$features'.")
  }

  df       <- umap_result$df
  feat_mat <- as.matrix(umap_result$features)
  storage.mode(feat_mat) <- "numeric"

  if (!all(c("V1", "V2") %in% colnames(df))) {
    stop("'umap_result$df' must contain columns 'V1' and 'V2'.")
  }
  if (!colour_by %in% colnames(df)) {
    stop("'colour_by' column '", colour_by, "' not found in umap_result$df.")
  }
  if (!sample_id_column %in% colnames(df)) {
    stop("'sample_id_column' ('", sample_id_column, "') not found in umap_result$df.")
  }

  # ── same filtering as make_umap_scatterplot ───────────────────────────────
  if (drop_composite) {
    df <- dplyr::filter(df, !is.na(.data[[colour_by]]),
                        !grepl("COMP", .data[[colour_by]]))
  }
  if (drop_other) {
    df <- dplyr::filter(df, !is.na(.data[[colour_by]]),
                        .data[[colour_by]] != "Other",
                        .data[[colour_by]] != "NOS")
  }
  if (high_confidence) {
    if (!"Confidence" %in% colnames(df)) {
      warning("'Confidence' column not found; ignoring high_confidence filter.")
    } else {
      df <- dplyr::filter(df, Confidence > 0.7)
    }
  }
  if (is.numeric(df[[colour_by]])) {
    df[[colour_by]] <- factor(df[[colour_by]])
  }

  if (!missing(custom_colours)) {
    cols <- custom_colours
  } else {
    cols <- GAMBLR.helpers::get_gambl_colours()
  }

  # ── build per-sample hover text ───────────────────────────────────────────
  hover_features <- vapply(df[[sample_id_column]], function(sid) {
    row    <- feat_mat[sid, , drop = TRUE]
    active <- names(row)[!is.na(row) & row > 0]
    if (length(active) == 0L) return("(no mutated features)")
    if (length(active) > max_features) {
      active <- c(active[seq_len(max_features)],
                  paste0("... +", length(active) - max_features, " more"))
    }
    paste(active, collapse = "<br>")
  }, character(1L))

  df$.hover_text <- paste0(
    "<b>", df[[sample_id_column]], "</b>",
    " (", df[[colour_by]], ")<br>",
    "<b>Mutated features:</b><br>",
    hover_features
  )

  # ── colour mapping ────────────────────────────────────────────────────────
  group_vals   <- unique(as.character(df[[colour_by]]))
  present_cols <- cols[intersect(names(cols), group_vals)]
  missing_grps <- setdiff(group_vals, names(present_cols))
  if (length(missing_grps) > 0) {
    message("No colour defined for group(s): ",
            paste(missing_grps, collapse = ", "),
            ". Assigning fallback colours.")
    fallback     <- grDevices::hcl.colors(length(missing_grps), palette = "Dark 2")
    present_cols <- c(present_cols, setNames(fallback, missing_grps))
  }

  # render "Other" beneath non-Other groups (plotly draws in factor order)
  level_order <- intersect(
    c("Other", setdiff(names(present_cols), "Other")),
    group_vals
  )
  df$.colour <- factor(df[[colour_by]], levels = level_order)
  df         <- dplyr::arrange(df, .colour)

  # ── build plotly chart ────────────────────────────────────────────────────
  plotly::plot_ly(
    data          = df,
    x             = ~V1,
    y             = ~V2,
    color         = ~.colour,
    colors        = present_cols[level_order],
    type          = "scatter",
    mode          = "markers",
    marker        = list(size = point_size, opacity = alpha),
    text          = ~.hover_text,
    hovertemplate = "%{text}<extra></extra>"
  ) %>%
    plotly::layout(
      title  = title,
      xaxis  = list(title = "UMAP 1", zeroline = FALSE, showgrid = FALSE),
      yaxis  = list(title = "UMAP 2", zeroline = FALSE, showgrid = FALSE),
      legend = list(title = list(text = paste0("<b>", colour_by, "</b>")))
    )
}
