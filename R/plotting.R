#' Heatmap visualization of mutations in nearest neighbors for a sample
#'
#' Generates a heatmap of feature values for the nearest neighbors
#' of a specified sample,
#' based on a DLBCLone model object. Supports outputs from:
#' - DLBCLone_KNN (with predict_unlabeled = TRUE)
#' - DLBCLone_optimize_params
#' - predict_single_sample_DLBCLone
#'
#' @param this_sample_id Character. The sample ID for which to plot
#' the nearest neighbor heatmap.
#' @param DLBCLone_model List. A DLBCLone model object as described above.
#' @param truth_column Character. Column name in predictions/metadata
#' to use as truth (default "lymphgen").
#' @param metadata_cols Optional character vector of additional metadata
#' columns to annotate on rows.
#' @param clustering_distance Distance for row clustering (default "binary").
#' @param font_size Numeric. Font size for labels (default 14).
#'
#' @return A ComplexHeatmap object (drawn).
#' @importFrom dplyr filter select left_join mutate pull
#' @importFrom tidyr separate pivot_longer
#' @importFrom tibble rownames_to_column column_to_rownames
#' @import ComplexHeatmap
#' @importFrom grid gpar
#' @importFrom circlize colorRamp2
#' @importFrom stats setNames
#' @export
nearest_neighbor_heatmap <- function(
  this_sample_id,
  DLBCLone_model,
  truth_column = "lymphgen",
  metadata_cols = NULL,
  clustering_distance = "binary",
  font_size = 14,
  gene_orientation = "column",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_legend = FALSE,
  draw = TRUE,
  keep_metafeatures = FALSE
){
  stopifnot(is.list(DLBCLone_model), "type" %in% names(DLBCLone_model))

  # Helper: get neighbor vector for a sample, handling both encodings
  .neighbors_for <- function(df_neighbors, sid) {
    if (!"sample_id" %in% names(df_neighbors)) return(character(0))
    row <- dplyr::filter(df_neighbors, .data$sample_id == sid)
    if (nrow(row) == 0) return(character(0))

    # Case A: already separated into N1..Nk columns
    n_cols <- grep("^N[0-9]+$", names(row), value = TRUE)
    if (length(n_cols) > 0) {
      vals <- unname(unlist(row[1, n_cols, drop = TRUE]))
      return(as.character(na.omit(vals)))
    }

    # Case B: single neighbor_id column "id1,id2,..."
    if ("neighbor_id" %in% names(row)) {
      vals <- strsplit(row$neighbor_id[1], ",", fixed = TRUE)[[1]]
      return(as.character(na.omit(vals)))
    }

    character(0)
  }

  # Decide model type + predicted label column + feature matrix source
  model_type <- DLBCLone_model$type
  pred_col   <- "DLBCLone_wo"
  greedy_col <- "DLBCLone_w"
  feats_mat  <- NULL
  neighbors  <- character(0)

  if (model_type %in% c("DLBCLone_KNN", "DLBCLone_predict")) {
    # Needs unlabeled_neighbors present (for KNN path)
    if (model_type == "DLBCLone_KNN") {
      if (!"unlabeled_neighbors" %in% names(DLBCLone_model) || is.null(DLBCLone_model$unlabeled_neighbors)) {
        message("No neighbors found (unlabeled_neighbors missing or NULL). Returning NULL.")
        return(NULL)
      }
      neighbors <- .neighbors_for(DLBCLone_model$unlabeled_neighbors, this_sample_id)
      pred_col  <- "DLBCLone_ko"
      greedy_col <- "DLBCLone_k"
    } else {
      # DLBCLone_predict stores neighbor_ids in prediction$neighbor_id (for non-self rows)
      # and unprocessed_votes$neighbor_id (for self rows)
      # Reconstruct a tiny neighbors df on the fly to reuse helper

      DLBCLone_model$predictions = DLBCLone_model$prediction

      pred_row <- DLBCLone_model$unprocessed_votes %>%
        dplyr::filter(.data$sample_id == .env$this_sample_id)
      if (nrow(pred_row) == 0 || !"neighbor_id" %in% names(pred_row)) {
        message("No neighbors found for sample ", this_sample_id, ". Returning NULL.")
        return(NULL)
      }
      neighbors <- strsplit(pred_row$neighbor_id[1], ",", fixed = TRUE)[[1]]
      neighbors <- as.character(na.omit(neighbors))
      pred_col <- "DLBCLone_wo" #eventually needs to support io or wo
      greedy_col <- "DLBCLone_w"
      pred_row <- DLBCLone_model$prediction %>%
        dplyr::filter(.data$sample_id == .env$this_sample_id)
      print(pred_row)
    }

    #feats_mat <- DLBCLone_model$features_df
    #if (is.null(feats_mat)) {
    #  stop("features_df is missing in the provided model object.")
    #}
    index_feats_mat <- DLBCLone_model$features #features for index sample
    other_feats_mat <- DLBCLone_model$umap_input #features for all training samples
    #check that index_feats_mat and other_feats_mat have the same columns
    if(!all(colnames(index_feats_mat) == colnames(other_feats_mat))){
      stop("features and umap_input must have the same columns")
    }

    feats_mat <- rbind(index_feats_mat, other_feats_mat)


  } else if (model_type == "DLBCLone_optimize_params") {
    # neighbors via predictions$neighbor_id for the sample
    if (!"predictions" %in% names(DLBCLone_model)) {
      stop("DLBCLone_optimize_params model missing $predictions")
    }

    pred_row <- DLBCLone_model$predictions %>%
      dplyr::filter(.data$sample_id == .env$this_sample_id)
    if (nrow(pred_row) == 0 || !"neighbor_id" %in% names(pred_row)) {
      message("No neighbors found for sample ", this_sample_id, ". Returning NULL.")
      return(NULL)
    }
    neighbors <- strsplit(pred_row$neighbor_id[1], ",", fixed = TRUE)[[1]]
    neighbors <- as.character(na.omit(neighbors))
    neighbor_dists <- strsplit(pred_row$distance[1], ",", fixed = TRUE)[[1]]
    neighbor_dists <- as.numeric(na.omit(neighbor_dists))
    other_neighbors <- strsplit(pred_row$other_neighbor[1], ",", fixed = TRUE)[[1]]
    other_neighbors <- as.character(na.omit(other_neighbors))
    neighbors = unique(c(neighbors, other_neighbors))
    feats_mat <- DLBCLone_model$features
    if (is.null(feats_mat)) {
      stop("features is missing in the optimize_params model object.")
    }



  } else {
    stop("DLBCLone_model$type must be one of: DLBCLone_KNN, DLBCLone_predict, DLBCLone_optimize_params")
  }

  # Make sure the focal sample itself is included (to show in heatmap/annotation)
  neighbors <- unique(c(this_sample_id, neighbors))
  if (length(neighbors) == 0) {
    message("No neighbors to display. Returning NULL.")
    return(NULL)
  }

  # Subset feature matrix; require that all neighbor rows exist

  xx <- feats_mat[neighbors, , drop = FALSE]

  if(!keep_metafeatures){
    xx <- xx %>% select(-ends_with("_feats"))
  }
  if (any(is.na(rownames(xx)))) {
    print(neighbors)
    stop("Some neighbor samples are missing from
      the feature matrix (features_df/features).")
  }
  # Keep only non-zero (informative) features
  xx <- xx[, colSums(xx) > 0, drop = FALSE]
  if (ncol(xx) == 0) {
    message("All selected features are zero across
      neighbors; nothing to plot. Returning NULL.")
    return(NULL)
  }

  # Colors for heatmap
  top <- max(xx, na.rm = TRUE)
  mid <- top / 2
  col_fun <- circlize::colorRamp2(c(0, mid, top), c("white", "#FFB3B3", "red"))

  # ----- Build row annotation data -----
  # Start from predictions; join metadata if present
  if(DLBCLone_model$type=="DLBCLone_predict"){
    preds <- DLBCLone_model$prediction
    
  }else{
    preds <- DLBCLone_model$predictions
  }

  if (is.null(preds) || !"sample_id" %in% names(preds)) {
    stop("Model object missing predictions with 'sample_id'.")
  }

  # If metadata present, join so we can pull truth_column from there if needed
  if ("metadata" %in% names(DLBCLone_model) && !is.null(DLBCLone_model$metadata)) {
    if(DLBCLone_model$type=="DLBCLone_predict"){
      #need to pool together training sample metadata with predictions
      train_meta = DLBCLone_model$metadata %>%
        dplyr::select(dplyr::all_of(c("sample_id", truth_column, metadata_cols))) %>%
        filter(!sample_id %in% preds$sample_id)
       #preds <- dplyr::left_join(preds, DLBCLone_model$metadata, by = "sample_id")

      preds <- dplyr::bind_rows(
        train_meta,
        preds
      ) %>% dplyr::distinct(sample_id, .keep_all = TRUE)

      #preds %>% select(sample_id,lymphgen,DLBCLone_wo) %>% head() %>% print()
    }else{
      preds <- dplyr::left_join(preds, DLBCLone_model$metadata, by = "sample_id")
    }
    #preds <- dplyr::left_join(preds, DLBCLone_model$metadata, by = "sample_id")
    #print(colnames(preds))
  }else{
    warning("No metadata found in model; truth_column and metadata_cols may be missing.")
  }

  # Columns to annotate
  cols_to_select <- unique(c("sample_id", truth_column, pred_col, greedy_col, metadata_cols))
  #print(cols_to_select)
  cols_to_select <- cols_to_select[cols_to_select %in% names(preds)]

  row_df <- preds %>%
    dplyr::select(dplyr::all_of(cols_to_select)) %>%
    dplyr::filter(.data$sample_id %in% rownames(xx))

  # If the focal sample is missing from predictions (common for "unlabeled"),
  # try unlabeled_predictions
  sample_class <- NA_character_
  if (!this_sample_id %in% row_df$sample_id) {
    if("prediction" %in% names(DLBCLone_model) &&
      DLBCLone_model$type == "predict_single_sample_DLBCLone"){
      stop(paste(this_sample_id,
        "not found in predictions; cannot look up class for predict_single_sample_DLBCLone"))
    }
    if ("unlabeled_predictions" %in% names(DLBCLone_model) &&
        !is.null(DLBCLone_model$unlabeled_predictions)) {
      extra <- DLBCLone_model$unlabeled_predictions %>%
        dplyr::select(dplyr::any_of(cols_to_select)) %>%
        dplyr::filter(.data$sample_id %in% rownames(xx))
      row_df <- dplyr::bind_rows(row_df, extra)

      if (this_sample_id %in% extra$sample_id && pred_col %in% names(extra)) {
        sample_class <- extra %>%
          dplyr::filter(.data$sample_id == .env$this_sample_id) %>%
          dplyr::pull(.data[[pred_col]]) %>% as.character()
      }
    } else {
      # Create a placeholder row with NAs for missing fields
      missing_row <- as.data.frame(setNames(rep(NA,
        length(cols_to_select)), cols_to_select))
      missing_row$sample_id <- this_sample_id
      row_df <- dplyr::bind_rows(row_df, missing_row)
    }
  } else {
    # Present in predictions; try to get predicted class cleanly
    if (pred_col %in% names(row_df)) {
      sample_class <- row_df %>%
        dplyr::filter(.data$sample_id == .env$this_sample_id) %>%
        dplyr::pull(.data[[pred_col]]) %>% as.character()
    }
    # For optimize_params, original code looked up
    # predicted_label_optimized explicitly
    if (model_type == "DLBCLone_optimize_params" &&
        pred_col %in% names(DLBCLone_model$predictions)) {
      sample_class <- DLBCLone_model$predictions %>%
        dplyr::filter(.data$sample_id == .env$this_sample_id) %>%
        dplyr::pull(.data[[pred_col]]) %>% as.character()
      sample_greedy_class <- DLBCLone_model$predictions %>%
        dplyr::filter(.data$sample_id == .env$this_sample_id) %>%
        dplyr::pull(.data[[greedy_col]]) %>% as.character()
    }
  }

  # Align rows to feature matrix order for annotation
  row_df <- row_df %>% tibble::column_to_rownames("sample_id")
  row_df <- row_df[rownames(xx), , drop = FALSE]

  # Build annotation colour map (reuse one palette; unknowns will recycle)
  anno_colours <- get_gambl_colours()
  anno_list <- list()
  if (truth_column %in% colnames(row_df)) anno_list[[truth_column]] <- anno_colours
  if (!is.null(pred_col) && pred_col %in% colnames(row_df)) anno_list[[pred_col]] <- anno_colours
  if (!is.null(greedy_col) && greedy_col %in% colnames(row_df)) anno_list[[greedy_col]] <- anno_colours
  if (!is.null(metadata_cols)) {
    for (mc in metadata_cols) {
      if (mc %in% colnames(row_df)) anno_list[[mc]] <- anno_colours
    }
  }
  # Title
  title_text <- paste("Sample", this_sample_id, "classified as",
    ifelse(is.na(sample_class), "<NA>", sample_class))

  if(gene_orientation == "column"){
    row_anno <- ComplexHeatmap::rowAnnotation(
        df  = row_df,
        col = anno_list,
        annotation_name_gp = grid::gpar(fontsize = font_size),
        show_legend = show_legend
      )
          # Draw heatmap
      ht <- ComplexHeatmap::Heatmap(
        xx,
        col = col_fun,
        right_annotation = row_anno,
        clustering_distance_rows = clustering_distance,
        column_title = title_text,

        column_names_gp = grid::gpar(fontsize = font_size),
        column_title_gp = grid::gpar(fontsize = font_size),
        row_names_gp = grid::gpar(fontsize = font_size),
        show_heatmap_legend = FALSE,
        cluster_rows = cluster_rows,
        cluster_columns = cluster_columns
      )

  }else{
    for_heatmap <- t(xx)
    annot_df <- row_df[colnames(for_heatmap), , drop = FALSE]

    # build a *column* annotation (same colours list you used for rows)
    col_anno <- ComplexHeatmap::HeatmapAnnotation(
      df  = annot_df,
      col = anno_list,  # keep if you want coloured discrete annotations
      annotation_name_gp = grid::gpar(fontsize = font_size),
      show_legend = show_legend,
      which = "column"
    )

    ht <- ComplexHeatmap::Heatmap(
      for_heatmap,
      col = col_fun,
      bottom_annotation = col_anno,                 # or top_annotation = col_anno
      clustering_distance_columns = clustering_distance,
      show_heatmap_legend = FALSE,
      column_title = title_text,
      column_title_gp = grid::gpar(fontsize = font_size),
      row_names_gp = grid::gpar(fontsize = font_size),
      column_names_gp = grid::gpar(fontsize = font_size),
      show_heatmap_legend=FALSE,
      cluster_rows = cluster_rows,
      cluster_columns = cluster_columns
    )

  }
  if(draw){
    ComplexHeatmap::draw(
      ht,
      heatmap_legend_side = "bottom",
      annotation_legend_side = "bottom"
    )
  }else{
    return(ht)
  }

}




#' Basic UMAP Scatterplot
#'
#' Generates a simple UMAP scatterplot for visualizing sample
#' clustering or separation.
#'
#' @param optimized Data frame containing at least V1, V2,
#' sample_id, and grouping columns.
#' @param plot_samples Optional character vector of sample_ids
#' to label in the plot
#' @param colour_by Column name to color points by. Defaults
#' to `truth_column`.
#' @param truth_column Name of the truth/ground-truth column
#' (default: "lymphgen").
#' @param pred_column  Name of the predicted-class column
#' (default: "DLBCLone_ko").
#' @param other_label  Label used for the outgroup/unclassified
#' class (default: "Other").
#' @param title Plot title.
#' @param use_plotly Logical; if FALSE and `plot_samples` provided,
#' draw static labels.
#' @param custom_colours Optional named vector of colors for groups;
#' falls back to `get_gambl_colours()`.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#'
#' \dontrun{
#' my_umap = make_and_annotate_umap(my_data, my_metadata)
#'
#' basic_umap_scatterplot(my_umap$df,
#'                        plot_samples = "some_sample_ID",
#'                        colour_by = "DLBCLone_ko")
#' }
basic_umap_scatterplot <- function(optimized,
                                   plot_samples = NULL,
                                   colour_by    = NULL,
                                   truth_column = "lymphgen",
                                   pred_column  = "DLBCLone_ko",
                                   other_label  = "Other",
                                   title        = "UMAP based on selected features",
                                   use_plotly   = TRUE,
                                   custom_colours = NULL) {
  stopifnot(all(c("V1","V2","sample_id") %in% colnames(optimized)))
  stopifnot(is.data.frame(optimized))
stopifnot(all(c("V1", "V2") %in% colnames(optimized)))
message("colour_by: ", colour_by)
  colour_by <- colour_by %||% truth_column

  # dynamic label text (e.g., "lymphgen" and "DLBCLone_ko")
  optimized_label <- optimized %>%
    mutate(
      label = paste0(
        sample_id,
        "\n", truth_column, ":\n", .data[[truth_column]],
        "  ", pred_column, ":\n", .data[[pred_column]]
      )
    )

  # points to label
  label_points <- dplyr::filter(optimized_label, sample_id %in% (plot_samples %||% character())) %>%
    mutate(label_x = V1 + 0.75, label_y = V2 - 0.75)

  # base mapping: color by chosen column
  aes_cols <- aes(
    x = V1, y = V2,
    sample_id = sample_id,
    color = .data[[colour_by]]
  )

  # draw outgroup first (by truth column) to keep it visually underneath
  p <- ggplot()
  p <- p +
    geom_point(
      data = dplyr::filter(optimized, .data[[truth_column]] == other_label),
      mapping = aes_cols
    ) +
    geom_point(
      data = dplyr::filter(optimized, .data[[truth_column]] != other_label),
      mapping = aes_cols
    )

  # palette
  pal <- custom_colours %||% get_gambl_colours()
  p <- p + scale_color_manual(values = pal) +
    guides(color = guide_legend(title = if (identical(colour_by, truth_column)) "Original Class" else "Predicted Class"))

  # optional static labels when not using plotly
  if (!is.null(plot_samples) && length(plot_samples) > 0 && !isTRUE(use_plotly)) {
    p <- p +
      geom_segment(
        data = label_points,
        aes(x = V1, y = V2, xend = label_x, yend = label_y),
        arrow = arrow(length = unit(0.02, "npc")),
        color = "black"
      ) +
      geom_label(
        data = label_points,
        aes(x = label_x, y = label_y, label = label),
        size = 3,
        fill = "white",
        label.size = 0.25
      )
  }

  p <- p +
    labs(title = title) +
    theme_minimal()

  return(p)
}


#' Summarize and Export DLBCLone Model Results
#'
#' Generates and saves a set of summary plots and tables for
#' a DLBCLone model, including UMAP scatterplots, alluvial
#' plots, and oncoplots.
#' Results are saved as PDF files in a directory named after the provided
#' base name.
#'
#' @param base_name Character. The base name (and directory)
#' for saving output files.
#' @param optimized_model List. The output from DLBCLone optimization,
#' containing predictions, features, and metadata.
#'
#' @details
#' - Creates a directory for results if it does not exist.
#' - Saves UMAP scatterplots for all samples and for non-"Other" samples.
#' - Generates alluvial plots for different DLBCLone predictions
#' - Exports an oncoplot summarizing mutation and classification results.
#' - Uses `make_umap_scatterplot`, `make_alluvial`, and `prettyOncoplot`
#' for visualization.
#'
#' @return No return value. Side effect: writes multiple PDF files to disk.
#'
#' @import ggalluvial grDevices
#'
#' @export
#'
#' @examples
#' \dontrun{
#' DLBCLone_summarize_model(optimized_model = all_features_optimized,
#'                          base_name="model_summaries/Lymphgen_all_features")
#' }
DLBCLone_summarize_model = function(base_name,
                                    optimized_model){
  base_dir = here::here()
  full_dir = paste0(base_dir,"/",base_name)
  if(!dir.exists(full_dir)){
    dir.create(full_dir)
  }
  #save the predictions
  pred_out = paste0(full_dir,"/DLBCLone_predictions_bulk.tsv")
  print(paste("Saving predictions to",pred_out))
  to_save = optimized_model$predictions %>%
    select(sample_id,!!sym(optimized_model$truth_column),DLBCLone_i:DLBCLone_wo,confidence:other_score,V1,V2)

  write_tsv(to_save,
                file = pred_out)
  umap1 = paste0(full_dir,"/UMAP_all.pdf")
  cairo_pdf(umap1,width=8,height=8)
  p = make_umap_scatterplot(optimized_model$df,colour_by = optimized_model$truth_column,
                            title="all samples, projected")
  print(p)
  dev.off()

  umap2 = paste0(full_dir,"/UMAP_no_Other.pdf")
  cairo_pdf(umap2,width=8,height=8)
  p = make_umap_scatterplot(optimized_model$df,
                            colour_by = optimized_model$truth_column,
                            drop_other = T,title="non-Other samples, projected")
  print(p)
  dev.off()
  cairo_pdf(paste0(full_dir,"/Alluvial_DLBCLone_io.pdf"),width=8,height=12)
  p = make_alluvial(optimized_model,pred_column = "DLBCLone_io",
                    truth_column = optimized_model$truth_column)
  print(p)
  dev.off()
  cairo_pdf(paste0(full_dir,"/Alluvial_DLBCLone_wo.pdf"),width=8,height=12)
  p = make_alluvial(optimized_model,pred_column = "DLBCLone_wo",
                    truth_column = optimized_model$truth_column)
  print(p)
  dev.off()
  cairo_pdf(paste0(full_dir,"/Oncoplot.pdf"),width=16,height=14)
  mc = c("lymphgen","DLBCLone_wo","DLBCLone_io","DLBCLone_w","DLBCLone_i")
  sc = c("DLBCLone_wo","DLBCLone_io","lymphgen","DLBCLone_w","DLBCLone_i","confidence")
  if("DLBClass" %in% colnames(optimized_model$predictions)){
    mc = c(mc,"DLBClass")
    sc = c(sc,"DLBClass")
  }
  mc = mc[mc %in% colnames(optimized_model$predictions)]
  sc = sc[sc %in% colnames(optimized_model$predictions)]
  prettyOncoplot(all_maf_with_s,
               these_samples_metadata =optimized_model$predictions,
               genes = colnames(optimized_model$features),
               minMutationPercent = 1,
               metadataColumns =  mc,
               sortByColumns = sc,
               numericMetadataColumns = "confidence",
               simplify_annotation = T,
               cluster_rows=T
               )
  dev.off()


}


#' @title Make Neighborhood Plot
#' @description
#' Generates a UMAP plot highlighting the neighborhood of a given sample,
#' showing its nearest neighbors and their group assignments.
#'
#' @param single_sample_prediction_output A list containing prediction
#' results and annotation data frames.
#'        Must include elements \code{prediction} (data frame with prediction
#' results) and \code{anno_df} (data frame with UMAP coordinates and
#' annotations).
#' @param this_sample_id Character. The sample ID for which the neighborhood
#' plot will be generated.
#' @param prediction_in_title Logical. If \code{TRUE}, includes the predicted
#' label in the plot title.
#' @param add_circle Plot will include a circle surrounding the set of
#' neighbors. Set to FALSE to disable.
#' @param label_column Specify the column that contains the DLBCLone
#' prediction you want to show.
#' Default: the optimized Other-balanced "DLBCLone_wo". 
#' @param point_size Numeric. The size of the points in the plot.
#' Default: 0.5.
#' @param point_alpha Numeric. The transparency level of the points
#' in the plot. Default: 0.9.
#' @param line_alpha Numeric. The transparency level of the lines
#' connecting the sample to its neighbors. Default: 0.9.
#' @return A \code{ggplot2} object representing the UMAP plot with
#' the selected sample and its neighbors highlighted.
#'
#' @details
#' The function extracts the nearest neighbors of the specified sample,
#' draws segments connecting the sample to its neighbors, and colors
#' points by group (e.g., lymphgen subtype). The plot title can optionally
#' include the predicted label.
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom rlang sym
#'
#' @export
#' @examples
#'
#' # Assuming 'pred_output' is the result of DLBCLone_predict_single_sample
#' # on sample_id "DLBCL10951T". 
#' \dontrun{
#'  make_neighborhood_plot(pred_output, "DLBCL10951T")
#' }
make_neighborhood_plot <- function(single_sample_prediction_output,
                                  this_sample_id,
                                  prediction_in_title = TRUE,
                                  add_circle = TRUE,
                                  label_column = "DLBCLone_wo",
                                  point_size = 0.5,
                                  point_alpha = 0.9,
                                  line_alpha = 0.9) {

  circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
    r = diameter / 2
    tt <- seq(0,2*pi,length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
  }
  
 
  if(missing(single_sample_prediction_output)){
    stop("single_sample_prediction_output is required")
  }
  truth_column = single_sample_prediction_output$truth_column
  if(missing(this_sample_id)){
    stop("this_sample_id is required")
  }
  if(!this_sample_id %in% single_sample_prediction_output$prediction$sample_id){
    stop(paste(this_sample_id,"not found in prediction data frame"))
  }
  training_predictions = single_sample_prediction_output$training_predictions 
  single_sample_prediction_output$prediction = filter(single_sample_prediction_output$prediction, sample_id==this_sample_id)

  xmin = min(training_predictions$V1, na.rm = TRUE)
  xmax = max(training_predictions$V1, na.rm = TRUE)
  ymin = min(training_predictions$V2, na.rm = TRUE)
  ymax = max(training_predictions$V2, na.rm = TRUE)
  #extract the sample_id for all the nearest neighbors with non-Other labels
  my_neighbours = filter(single_sample_prediction_output$prediction,
                         sample_id == this_sample_id) %>%
                  pull(neighbor_id) %>% strsplit(.,",") %>% unlist()

  #set up links connecting each neighbor to the sample's point
  links_df = filter(training_predictions,sample_id %in% my_neighbours) %>% 
    mutate(group=!!sym(truth_column))
  my_x = filter(single_sample_prediction_output$projection,
                sample_id==this_sample_id) %>% pull(V1)
  my_y = filter(single_sample_prediction_output$projection,
                sample_id==this_sample_id) %>% pull(V2)
  if(prediction_in_title){
    title = paste(this_sample_id,
                  pull(single_sample_prediction_output$prediction,
                       !!sym(label_column)))
  }else{
    title = this_sample_id
  }
  links_df = mutate(links_df,my_x=my_x,my_y=my_y)
  links_df = links_df %>% select(V1,V2,my_x,my_y,group) %>%
    mutate(length = abs(V1-my_x)+abs(V2-my_y))


  pp=ggplot(mutate(training_predictions,group=!!sym(truth_column)),
         aes(x=V1,y=V2,colour=group)) +
    geom_point(alpha=point_alpha,size=point_size) +
    geom_segment(data=links_df,aes(x=V1,y=V2,xend=my_x,yend=my_y),alpha=line_alpha)+
    scale_colour_manual(values=get_gambl_colours()) +
    ggtitle(title)+
    xlim(c(xmin,xmax)) +
    ylim(c(ymin,ymax)) +
    theme_minimal()
  if(add_circle){
    #add a circle around the sample
    max_d = 1
    d = max(links_df$length)*2
    if(d>max_d){
      d = max_d
    }
    circle = circleFun(c(my_x,my_y),diameter=d,npoints=100)
    pp = pp + geom_path(data=circle,aes(x=x,y=y),colour="black",alpha=1,inherit.aes=FALSE)
  }
  return(pp)
}


#' Make UMAP scatterplot
#'
#' @param df Data frame containing UMAP coordinates (V1, V2) and grouping columns.
#' Typically this is the `df` element from the output of `make_and_annotate_umap()`.
#' @param drop_composite Logical; if TRUE, drops samples with
#' composite labels (e.g., "EZB-COMP").
#' @param colour_by Column name to color points by (default: "lymphgen").
#' @param drop_other Logical; if TRUE, drops (hides) samples with
#' "Other" labels.
#' @param high_confidence Logical; if TRUE, only includes samples
#' with high confidence predictions.
#' @param custom_colours Named vector of custom colors for each group.
#' @param add_labels Logical; if TRUE, adds labels to points.
#' @param title Plot title.
#' @param base_size Base font and point size for the plot
#' (passed to theme_Morons()).
#' @param alpha Point transparency (default: 0.8).
#'
#' @returns A ggplot object.
#'
#' @import ggside
#' @export
#'
#' @examples
#' 
#' \dontrun{
#' my_umap = make_and_annotate_umap(my_data, my_metadata)
#' p = make_umap_scatterplot(my_umap$df,
#'                           colour_by = "lymphgen",
#'                           drop_other = TRUE)
#' print(p)
#' }
make_umap_scatterplot = function(df,
                                 drop_composite = TRUE,
                                 colour_by="lymphgen",
                                 drop_other = FALSE,
                                 high_confidence = FALSE,
                                 custom_colours,
                                 add_labels = FALSE,
                                 title = NULL,
                                 base_size=8,
                                 alpha = 0.8){
  xmin = min(df$V1,na.rm=TRUE)
  xmax = max(df$V1,na.rm=TRUE)
  ymin = min(df$V2,na.rm=TRUE)
  ymax = max(df$V2,na.rm=TRUE)
  if(!"cohort" %in% colnames(df)){
    df$cohort = "Unknown"
  }
  if(!missing(custom_colours)){
    cols = custom_colours
  }else{
    cols = get_gambl_colours()
  }
  if(drop_composite){
    df = filter(df,!is.na(!!sym(colour_by)),!grepl("COMP",!!sym(colour_by)))
  }
  if(drop_other){
    df = filter(df,!is.na(!!sym(colour_by)),!!sym(colour_by) != "Other",!!sym(colour_by) !="NOS")
  }
  if(high_confidence){
    df = filter(df,Confidence > 0.7)
  }
  if(is.numeric(df[[colour_by]])){
    df[[colour_by]]=factor(df[[colour_by]])
  }
  if(add_labels){
    labels = group_by(df,!!sym(colour_by)) %>%
      summarise(median_x = median(V1),median_y = median(V2))
  }
  unique_lg = unique(df[[colour_by]])

  if(any(!unique_lg %in% names(cols))){
    missing = unique_lg[!unique_lg %in% names(cols)]
    print(paste("missing colour for:",paste(missing,collapse=",")))

  }
  point_size = base_size/12
  p = ggplot(df,
             aes(x=V1,y=V2,colour=!!sym(colour_by),label=cohort)) +
             geom_point(alpha=0,size=point_size) +
             geom_point(data=df %>% filter(!!sym(colour_by) =="Other"),alpha=alpha,size=point_size) +
             geom_point(data=df %>% filter(!!sym(colour_by) !="Other"),alpha=alpha,size=point_size) +
    scale_colour_manual(values=cols) +
    scale_fill_manual(values=cols) +
    theme_Morons(base_size = base_size) +
    guides(colour = guide_legend(nrow = 1)) +
    xlim(xmin,xmax) +
    ylim(ymin,ymax) +
    theme(axis.title.x = element_blank(),
      axis.title.y = element_blank())
  if(!is.null(title)){
    p = p + ggtitle(title)
  }
  if(add_labels){
    p = p + geom_label_repel(data=labels,aes(x=median_x,y=median_y,label=!!sym(colour_by)))
  }
  #p = ggMarginal(p,groupColour = TRUE,groupFill=TRUE)
  p <- p +

  ggside::geom_xsidedensity(aes(y = after_stat(density), fill = !!sym(colour_by)),
                            alpha = 0.3, size = 0.2) +
  ggside::geom_ysidedensity(aes(x = after_stat(density), fill = !!sym(colour_by)),
                            alpha = 0.3, size = 0.2) +
  ggside::scale_xsidey_continuous(breaks = NULL, minor_breaks = NULL) +
  ggside::scale_ysidex_continuous(breaks = NULL, minor_breaks = NULL)
  #theme(axis.ticks.x = NULL)
  return(p)
}



#' Calculate Classification Accuracy and Per-Class Metrics based on Predictions
#'
#' Computes overall accuracy, balanced accuracy, and sensitivity
#' for predicted vs. true class labels.
#' Optionally excludes samples assigned to the "Other" class from
#' accuracy calculations.
#'
#' @param predictions Data frame containing predicted and true class labels.
#' @param truth Name of the column with true class labels (default: "lymphgen").
#' @param pred Name of the column with predicted class labels (default: "predicted_label").
#' @param per_group Logical; if TRUE, computes per-group accuracy metrics.
#' @param metric Character; type of accuracy to report ("accuracy" supported).
#'
#' @return A list with:
#'   \item{no_other}{Accuracy excluding samples assigned to "Other"}
#'   \item{per_class}{Average of balanced accuracy values for each class}
#'   \item{per_class_sensitivity}{Sensitivity per class}
#'   \item{overall}{Overall accuracy including all samples}
#'
#' @details
#' - Uses confusion matrices to compute accuracy metrics.
#' - Excludes "Other" class for no_other accuracy.
#' - Returns per-class metrics for further analysis.
#'
#' @examples
#' \dontrun{
#' result <- report_accuracy(predictions_df)
#' result$overall
#' result$per_class
#' }
#' @export
report_accuracy <- function(predictions,
                            truth = "lymphgen",
                            pred = "DLBCLone_io",
                            per_group = FALSE,
                            metric = "accuracy",
                            verbose = FALSE,
                            drop_other = TRUE,
                            skip_F1 = FALSE) {
  #Helper function to compute macro F1 score
  macro_f1 <- function(truth, pred, drop_other = drop_other, other_label = "Other", na_rm = TRUE) {
    stopifnot(length(truth) == length(pred))
    truth <- factor(truth)
    pred  <- factor(pred, levels = levels(truth))  # align levels

    classes <- levels(truth)
    if (drop_other && other_label %in% classes) {
      classes <- setdiff(classes, other_label)
    }

    f1_per_class <- sapply(classes, function(cls) {
      tp <- sum(truth == cls & pred == cls)
      fp <- sum(truth != cls & pred == cls)
      fn <- sum(truth == cls & pred != cls)

      prec <- if ((tp + fp) == 0) NA_real_ else tp / (tp + fp)
      rec  <- if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)

      if (is.na(prec) || is.na(rec) || (prec + rec) == 0) {
        if (na_rm) return(NA_real_) else return(0)  # choose policy
      }
      2 * prec * rec / (prec + rec)
    })

   mean(f1_per_class, na.rm = na_rm)

  }
  all_classes <- unique(predictions[[truth]])
  no_other_true <- filter(predictions, !!sym(truth) != "Other")
  no_other_pred <- filter(predictions, !!sym(pred) != "Other")

  classification_rate = nrow(no_other_pred) / nrow(predictions)
  #print(paste(classification_rate,"=", nrow(no_other_pred) ,"/" ,nrow(predictions)))
  if (verbose) {
    print(paste(nrow(no_other_true), "non-Other samples were assigned a label"))
    print(paste(
      filter(no_other_pred, !!sym(pred) != "Other") %>% nrow(),
      "non-Other samples were assigned to a non-Other class"))
  }


  conf_matrix_no <- confusionMatrix(
    factor(no_other_true[[pred]], levels = all_classes),
    factor(no_other_true[[truth]], levels = all_classes)
  )

  overall <- conf_matrix_no$overall[["Accuracy"]]

  if (metric == "accuracy") {
    acc_no <- overall
  } else {
    stop("unsupported metric")
  }

  conf_matrix <- confusionMatrix(
    factor(predictions[[pred]], levels = unique(predictions[[truth]])),
    factor(predictions[[truth]], levels = unique(predictions[[truth]]))
  )
  if (metric == "accuracy") {
    bal_acc <- conf_matrix$byClass[, "Balanced Accuracy"]
    sensitivity <- conf_matrix$byClass[, "Sensitivity"]
  } else {
    stop("unsupported metric")
  }
  mba <- mean(bal_acc, na.rm = TRUE)
  if(!skip_F1){
    f1 = macro_f1(predictions[[truth]], predictions[[pred]], drop_other = TRUE)
    harm = 4 / ( 1/acc_no + 1/f1 + 1/mba + 1/classification_rate)
  }else{
    f1 = NA
    harm = NA
  }


  overall <- conf_matrix$overall[["Accuracy"]]
  return(list(
    no_other = acc_no,
    accuracy_no_other = acc_no,
    macro_f1 = f1,
    harmonic_mean = harm,
    per_class = bal_acc,
    mean_balanced_accuracy = mba,
    classification_rate = classification_rate,
    per_class_sensitivity = sensitivity,
    overall = overall,
    confusion_matrix_no_other = conf_matrix_no,
    confusion_matrix = conf_matrix,
    num_unclass = nrow(no_other_pred)
  ))
}


#' Create an Alluvial Plot Comparing Original and Predicted Classifications
#'
#' This function generates a detailed alluvial plot to visualize the concordance
#' and discordance between original (e.g., Lymphgen) and predicted
#' (e.g., DLBCLone) class assignments for samples.
#' It supports annotation of concordance rates, per-group accuracy,
#' unclassified rates, and flexible labeling and coloring options.
#'
#' @param optimized List containing prediction results and metadata,
#' typically output from a DLBCLone optimization function.
#' @param pred Name of the column in predictions to use for the predicted
#' class (default: "predicted_label_optimized").
#' @param count_excluded_as_other Logical; if TRUE, samples excluded due
#' to missing features are counted as "Other" in the plot.
#' @param title Plot title (default: empty string).
#' @param group_order Character vector specifying the order of groups/classes
#' for axes and coloring.
#' @param add_accuracy_to_title Logical; if TRUE, adds accuracy/concordance
#' rate to the plot title.
#' @param accuracy_per_group Logical; if TRUE, computes and displays per-group accuracy.
#' @param accuracy_type Type of accuracy to report (default: "sensitivity").
#' @param original_name Name for the original class column (default: "Lymphgen").
#' @param new_name Name for the predicted class column (default: "DLBCLone").
#' @param nudge Amount to nudge labels horizontally (default: 0.03).
#' @param box_nudge Amount to nudge label boxes (default: 0.15).
#' @param min_flip_n Minimum number of samples for a flow to be labeled (default: 10).
#' @param add_percent Logical; if TRUE, adds percent concordance to labels.
#' @param add_unclass_rate Logical; if TRUE, adds unclassified rate annotation to the plot.
#' @param denom Denominator for stratum/flow width scaling (default: 20).
#' @param label_size Size of label text (default: 3).
#' @param label_lock_y Logical; if TRUE, locks label movement to the x-axis only.
#' @param label_line_flip_colour Logical; controls color assignment for label lines.
#' @param label_box_flip_colour Logical; controls color assignment for label boxes.
#' @param concordant_label_relative_pos Position for concordant labels:
#' 0 (left), 0.5 (middle), or 1 (right).
#' @param rotate
#' 
#' @return A ggplot2 object representing the alluvial plot.
#'
#' @details
#' - Visualizes flows between original and predicted classes,
#' highlighting concordant and discordant assignments.
#' - Annotates concordance rate, per-group accuracy, and unclassified
#' rate as specified.
#' - Supports flexible labeling, coloring, and axis ordering for
#' publication-quality plots.
#'
#' @examples
#' \dontrun{
#'   make_alluvial(optimized_result)
#' }
#' @import ggrepel
#' @export
make_alluvial <- function(
    optimized,
    count_excluded_as_other = FALSE,
    title = "",
    group_order = NULL,
    add_accuracy_to_title = TRUE,
    accuracy_per_group = TRUE,
    accuracy_type = "sensitivity",
    truth_column = "lymphgen",
    pred_name = "DLBCLone",
    truth_name = "Lymphgen",
    pred_column = "predicted_label_optimized",
    nudge = 0.03,
    box_nudge = 0.15,
    min_flip_n = 10,
    add_percent = TRUE,
    add_unclass_rate = TRUE,
    denom = 20,
    label_size=3,
    label_lock_y=TRUE,
    label_line_flip_colour=TRUE,
    label_box_flip_colour=TRUE,
    concordant_label_relative_pos = 0.5, #we only handle 0, 0.5 and 1
    verbose = FALSE,
    custom_colours = NULL,
    hide_legend = TRUE,
    rotate = FALSE
) {
  predictions = optimized$predictions
  if(!truth_column %in% colnames(predictions)){
    stop(paste("truth_column",truth_column,"not found in predictions"))
  }
  if(!pred_column %in% colnames(predictions)){
    stop(paste("pred_column",pred_column,"not found in predictions"))
  }

  if(is.null(group_order)){
    if("truth_classes" %in% names(optimized)){
      group_order = optimized$truth_classes
    }else{
      group_order = sort(unique(c(predictions[[pred_column]],predictions[[truth_column]])))
    }

    #print("setting group order:")
    #print(group_order)
  }




  if (accuracy_per_group) {
    accuracies <- report_accuracy(predictions,
        truth = truth_column,
        pred = pred_column,
        per_group = accuracy_per_group,
        skip_F1 = TRUE)
  }else{
    add_accuracy_to_title = FALSE
    add_percent = FALSE
  }
  if(count_excluded_as_other){
    excluded_meta = optimized$sample_metadata_no_features %>%
      mutate(!!pred_column := "Other")
    predictions = bind_rows(excluded_meta,predictions)
    full_denominator = nrow(predictions)

  } else if("total_samples_available" %in% names(optimized)){
      full_denominator = optimized$total_samples_available
  }else{
    full_denominator = nrow(predictions)
  }

  xx <- predictions %>%
    rename(
      !!truth_name := !!sym(truth_column)

    )

  xx <- xx %>% rename (
    !!pred_name := !!sym(pred_column))

  xx <- xx %>%
    group_by(!!sym(truth_name), !!sym(pred_name)) %>%
    summarize(num = n(), .groups = "drop")


  grid <- expand_grid(
    !!sym(truth_name) := group_order,
    !!sym(pred_name) := group_order
  )

  xx <- grid %>%
    left_join(xx, by = c(truth_name, pred_name)) %>%
    mutate(num = replace_na(num, 0))

  xx[[truth_name]] <- factor(xx[[truth_name]], levels = group_order)
  xx[[pred_name]] <- factor(xx[[pred_name]], levels = group_order)

  xx = xx %>% mutate(pair = paste(!!sym(truth_name),!!sym(pred_name),sep="-"))
  xy = filter(xx,!!sym(truth_name)!="Other") %>%
    mutate(match=ifelse(!!sym(truth_name)==!!sym(pred_name),TRUE,FALSE)) %>%
    group_by(match) %>%
    summarise(concordant=sum(num)) %>%
    ungroup() %>%
    mutate(total=sum(concordant)) %>%
    mutate(percent=100*concordant/total)

  pc_conc = filter(xy,match==TRUE) %>% pull(percent) %>% round(.,1)
  if("Other" %in% group_order){
    if(add_accuracy_to_title){
      title = paste0(title," Concordance (non-Other): ",pc_conc,"%")
    }
  }
  if(accuracy_per_group){
    bacc = round(accuracies$mean_balanced_accuracy,4)
  }
  

  # One side sort

  xx <- xx %>%
    group_by(!!sym(truth_name)) %>%
    arrange(!!sym(pred_name), .by_group = TRUE) %>%
    ungroup() %>%
    mutate(flow_id = row_number(), flow_label = ifelse(num > min_flip_n, as.character(num), ""))

  xx = xx %>%
    group_by(!!sym(pred_name)) %>%
    arrange(!!sym(truth_name),.by_group=TRUE) %>%
    ungroup() %>%
    mutate(flow_id2 = row_number())

  lodes <- ggalluvial::to_lodes_form(
    xx,
    key = "axis",
    axes = c(truth_name, pred_name),
    id = "flow_id"
  ) %>%
    mutate(
      axis = as.numeric(axis),
      nudge_dir = 0.1
    ) %>%
    left_join(xx %>% select(flow_id, !!sym(truth_name), !!sym(pred_name)), by = "flow_id")


  stratum_bases <- lodes %>%
    mutate(stratum = factor(stratum, levels = group_order)) %>%
    group_by(axis, stratum) %>%
    summarize(total = sum(num), .groups = "drop") %>%
    arrange(axis, desc(stratum)) %>%
    group_by(axis) %>%
    mutate(stratum_base = cumsum(total) - total) %>%
    ungroup()

  left_lodes <- lodes %>%
    left_join(stratum_bases, by = c("axis", "stratum")) %>%
    separate(pair,into=c("left_group","right_group"),sep = "-") %>%
    mutate(right_group = factor(right_group,levels=levels(stratum))) %>%
    mutate(left_group = factor(left_group,levels=levels(stratum))) %>%
    group_by(axis, stratum) %>%
    arrange(
      desc(stratum),
      desc(right_group),
      num,
      .by_group = TRUE) %>%
    mutate(
      flow_bottom = cumsum(num) - num,
      flow_y = stratum_base + flow_bottom + num / 2
    ) %>%
    ungroup()




  right_lodes <- lodes %>%
    left_join(stratum_bases, by = c("axis", "stratum")) %>%
    separate(pair,into=c("left_group","right_group"),sep = "-") %>%
    mutate(right_group = factor(right_group,levels=levels(stratum))) %>%
    mutate(left_group = factor(left_group,levels=levels(stratum))) %>%
    group_by(axis, stratum) %>%
    arrange(

      desc(right_group),
      desc(left_group),
      num,
      .by_group = TRUE) %>%
    mutate(
      flow_bottom = cumsum(num) - num,
      flow_y = stratum_base + flow_bottom + num / 2
    ) %>%
    ungroup()



  lodes_mean_left = group_by(left_lodes,right_group,left_group) %>%
    mutate(y=mean(flow_y)) %>% ungroup() %>%
    arrange(right_group,left_group) %>%
    mutate(nudge_dir = -nudge)

  lodes_mean_right = group_by(right_lodes,left_group,right_group) %>%
    mutate(y=mean(flow_y)) %>% ungroup() %>%
    arrange(left_group,right_group) %>%
    mutate(nudge_dir = nudge)

  if(concordant_label_relative_pos == 0){
    lodes_mean = lodes_mean_left %>% mutate(y=flow_y) %>%
      filter(axis==1)
    #position on far left
  }else if(concordant_label_relative_pos == 1){

    lodes_mean = lodes_mean_right  %>%
      mutate(y=flow_y) %>%
      filter(axis==2)
    #position on far right
  }else if(concordant_label_relative_pos == 0.5){
    lodes_mean = bind_rows(filter(right_lodes,axis==2),filter(left_lodes,axis==1)) %>%
      group_by(left_group,right_group) %>%
      mutate(y=mean(flow_y)) %>% ungroup() %>%
      arrange(left_group,right_group) %>%
      filter(axis==2) %>%
      mutate(nudge_dir = 0)
    #middle

  }else{
    stop("unhandled value for concordant_label_relative_pos. Provide 0 (left), 0.5 (middle) or 1 (right)")
  }

lodes_right  = filter(lodes_mean_right,axis==2,left_group!=right_group) %>%
  arrange(flow_id2) %>%
  mutate(nudge_dir = -nudge) %>%
  mutate(left_char=as.character(left_group),right_char = as.character(right_group))


lodes_left = filter(lodes_mean_left,axis==1,left_group!=right_group) %>%
  arrange(flow_id) %>%
  mutate(left_char=as.character(left_group),right_char = as.character(right_group))  %>%
  ungroup()

lodes_left = lodes_left %>%
  mutate(segment_colour = if (!label_line_flip_colour) left_char else right_char) %>%
  mutate(label_fill = if(label_box_flip_colour) right_group else left_group) %>%
  mutate(nudge_dir = nudge)

lodes_right = lodes_right %>%
  mutate(segment_colour = if (label_line_flip_colour) left_char else right_char) %>%
  mutate(label_fill = if(!label_box_flip_colour) right_group else left_group) %>%
  mutate(nudge_dir = -nudge)

lodes_match = filter(lodes_mean,left_group==right_group) %>%
  mutate(nudge_dir=0)

xx_denom <- predictions %>%
  group_by(!!sym(truth_column)) %>%
  summarize(denom = n(), .groups = "drop") %>%
  rename(
    stratum = !!sym(truth_column)
  )
if(add_percent){
    lodes_match = left_join(lodes_match,xx_denom,by="stratum") %>%
      mutate(percent = round(100*num/denom,1)) %>%
      filter(!is.na(denom)) %>%
      mutate(flow_label = paste0(left_group, ": ", num," (",percent,"%",")"))

  }
  if(add_unclass_rate){
    right_other = filter(lodes_mean_right,right_group=="Other",axis==2) %>%
      summarize(total=sum(num)) %>%
      mutate(unclass= round(100 * total / full_denominator,1)) %>%
      mutate(label=paste0(total,"\n (",unclass,"%)")) %>%
               mutate(!!sym(truth_name) := "Other",!!sym(pred_name) := "Other")
    unclass = pull(right_other,unclass)
    if(add_accuracy_to_title){
      title = paste0(title,", Classification rate: ",100-unclass,"%")
      title = paste0(title," Bal Acc: ",bacc)
    }
    left_other = filter(lodes_mean_left,left_group=="Other",axis==2) %>%
      summarize(total=sum(num)) %>%
      mutate(unclass= round(100 * total / full_denominator,1)) %>%
      mutate(label=paste0(total,"\n (",unclass,"%)")) %>%
      mutate(!!sym(truth_name) := "Other",!!sym(pred_name) := "Other")
  }

  if(!is.null(custom_colours)){
    custom_colours = custom_colours[unique(names(custom_colours))]
    custom_colours = custom_colours[names(custom_colours) %in% group_order]
    if(length(custom_colours) < length(group_order)){
      missing = group_order[!group_order %in% names(custom_colours)]
      message(paste("Missing colours for:",paste(missing,collapse=",")))
      print(names(custom_colours))
      stop("No custom colour provided for some of groups in group_order")
    }
    group_colours = custom_colours

  }else{
    group_colours = get_gambl_colours()
  }
  pp = ggplot(xx, aes(
    axis1 = !!sym(truth_name),
    axis2 = !!sym(pred_name),
    y = num
  )) +
    geom_alluvium(aes(fill = !!sym(truth_name)), width = 1 / denom) +
    geom_stratum(width = 1 / denom, aes(fill = after_stat(stratum)), color = "black") +
    geom_label(
      data = filter(lodes_match,left_group==right_group,flow_label!=""),
      stat = "identity",
      inherit.aes = FALSE,
      colour="white",
      size=label_size,
      aes(
        x = concordant_label_relative_pos + 1,
        y,
        label = flow_label,
        fill = right_group
      )

    ) 
    if(!rotate){
    pp = pp + 
    geom_label_repel(
      data = filter(lodes_right,flow_label!=""),
      min.segment.length = 0,
      stat = "identity",
      inherit.aes = FALSE,
      direction = ifelse(label_lock_y,"x","both"),
      colour="white",
      size=label_size,
      nudge_x = box_nudge,

      aes(
        x = axis + nudge_dir,
        flow_y,
        label = flow_label,
        fill = label_fill,
        segment.colour = segment_colour
      )

    ) +
    geom_label_repel(
      data = filter(lodes_left,flow_label!=""),
      min.segment.length = 0,
      stat = "identity",
      inherit.aes = FALSE,
      direction = ifelse(label_lock_y,"x","both"),
      colour="white",
      size=label_size,
      nudge_x = -box_nudge,
      aes(
        x = axis + nudge_dir,
        flow_y,
        label = flow_label,
        fill = label_fill,

        segment.colour = segment_colour

      )

    )
    }
    pp = pp +
    scale_x_discrete(limits = c(truth_name, pred_name), expand = c(.1, .05)) +
    scale_fill_manual(values = group_colours) +
    scale_colour_manual(values = group_colours,
                        aesthetics = c("color", "segment.color")) +

    #guides(color = "none",segment_color = "none") +
    labs(
      y = "Number of Samples",
      x = NULL
    ) +
    theme_minimal() +
    ggtitle(title)

  if(rotate){
    pp = pp + theme_minimal() + coord_flip()
  }
  if(hide_legend){
    pp = pp +
      theme(legend.position = "none")

  }
  if(add_unclass_rate){
    pp = pp +
      geom_label_repel(data=left_other,
                       direction = "x",
                       x=1,y=20,
                 aes(fill=!!sym(truth_name),label=label),size=label_size) +
      geom_label_repel(data=right_other,
                       direction = "x",
                       x=2,y=20,
                 aes(fill=!!sym(pred_name),label=label),size=label_size)
  }
  pp + guides(
    fill = guide_legend(override.aes = list(label = "", colour = NA)),
  color = "none",
  segment.colour = "none",
  fill = guide_legend(title = "Class"))
}

#' Plot the result of a DLBCLone classification
#'
#' @param test_df Data frame containing the test data with
#' UMAP coordinates
#' @param train_df Data frame containing the training data
#' with UMAP coordinates
#' @param predictions_df Data frame containing the predictions
#' with UMAP coordinates
#' @param other_df Data frame containing the predictions for
#' samples in the "Other" class
#' @param details Single-row data frame with the best parameters
#' from DLBCLone_optimize_params
#' @param annotate_accuracy Set to true to add labels with
#' accuracy values
#' @param classes Vector of classes that were used in the
#' training and testing
#' @param label_offset Length of the label offset for the
#' accuracy labels
#' @param title1 additional argument
#' @param title2 additional argument
#' @param title3 additional argument
#'
#' @returns a ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' #add the dataset name to the metadata if it's not already
#' #there (required for the plot work)
#' lymphgen_A53_DLBCLone$df$dataset = "GAMBL"
#'
#' DLBCLone_train_test_plot(
#'  test_df = lymphgen_A53_DLBCLone$df,
#'  train_df = lymphgen_A53_DLBCLone$df,
#'  predictions_df = lymphgen_A53_DLBCLone$predictions,
#'  #other_df = lymphgen_A53_DLBCLone$predictions_other,
#'  #required only when "Other" was in the truth_classes
#'  details = lymphgen_A53_DLBCLone$best_params,
#'  classes = c("MCD","EZB","BN2","ST2","N1","A53","Other"),
#'  annotate_accuracy=TRUE,label_offset = 1)
#' }
DLBCLone_train_test_plot = function(test_df,
                           train_df,
                           predictions_df,
                           other_df,
                           details,
                           annotate_accuracy = FALSE,

                           classes = c("BN2","ST2","MCD","EZB","N1"),
                           label_offset = 2,
                           title1="Original Class",
                           title2="DLBCLone Predicted Class",
                           title3 ="DLBCLone Predicted Class (Other)",base_size = 1){
  title = ""
  if(!missing(details)){
    title = paste0("N_class:",details$num_classes," N_feats:",details$num_features," k=",details$k," threshold=",details$threshold," bacc=",round(details$accuracy,3))

  }
  if(annotate_accuracy){
    if("BN2" %in% classes){
      acc_df = data.frame(lymphgen = classes,
                          accuracy = c(
                            details$BN2_bacc,
                            details$EZB_bacc,
                            details$MCD_bacc,
                            details$ST2_bacc,
                            details$Other_bacc,
                            details$A53_bacc))
    }else if("C1" %in% classes){
      acc_df = data.frame(lymphgen = c("C1","C2","C3","C4","C5"),
                          accuracy = c(details$C1_bacc,
                                       details$C2_bacc,
                                       details$C3_bacc,
                                       details$C4_bacc,
                                       details$C5_bacc))
    }else{
      stop("no labels to add?")
    }

  }
  # Add the predicted labels for Other (unclassified) cases, if provided
  if(!missing(other_df)){
    in_df = bind_rows(train_df,
                      test_df,
                      mutate(predictions_df,dataset=title2,lymphgen=predicted_label),
                      mutate(other_df,dataset=title3,lymphgen=predicted_label)
                      )
    in_df = mutate(in_df,dataset = factor(dataset,levels=unique(c(unique(train_df$dataset),title1,title2,title3))))
  }else{
    in_df = bind_rows(train_df,
                      test_df,
                      mutate(predictions_df,dataset=title2,lymphgen=predicted_label)
                      )
    in_df = mutate(in_df,dataset = factor(dataset,levels=unique(c(unique(train_df$dataset),title1,title2))))
  }

  pp = ggplot(in_df) +
    geom_point(aes(x=V1,y=V2,colour=lymphgen),alpha=0.8) +
    scale_colour_manual(values=get_gambl_colours()) +
    facet_wrap(~dataset,ncol=1) +
    theme_Morons(base_size=base_size) + ggtitle(title)
  if(annotate_accuracy){
    #add labels and set nudge direction based on what quadrant each group sits in
    centroids = filter(predictions_df,predicted_label %in% classes) %>%
      group_by(predicted_label) %>%
      summarise(mean_V1=median(V1),mean_V2=median(V2)) %>%
      mutate(nudge_x=sign(mean_V1),nudge_y = sign(mean_V2)) %>%
      mutate(lymphgen=predicted_label)
    #print(centroids)
    centroids = left_join(centroids,acc_df) %>%
      mutate(label=paste(lymphgen,":",round(accuracy,3)))

    centroids$dataset = title2
    pp = pp + geom_label_repel(data=filter(centroids,nudge_y < 0, nudge_x < 0),
                              aes(x=mean_V1,y=mean_V2,label=label),fill="white",size=5,nudge_y = -1 * label_offset , nudge_x = -1 * label_offset) +
      geom_label_repel(data=filter(centroids,nudge_y < 0, nudge_x > 0),
                      aes(x=mean_V1,y=mean_V2,label=label),size=5,nudge_y = -1 * label_offset , nudge_x = 1 * label_offset) +
      geom_label_repel(data=filter(centroids,nudge_y > 0, nudge_x < 0),
                      aes(x=mean_V1,y=mean_V2,label=label),size=5,nudge_y = 1 * label_offset , nudge_x = -1 * label_offset) +
      geom_label_repel(data=filter(centroids,nudge_y > 0, nudge_x > 0),
                      aes(x=mean_V1,y=mean_V2,label=label),fill="white",size=5,nudge_y = 1 * label_offset , nudge_x = 1 * label_offset)
  }
  pp + guides(colour = guide_legend(nrow = 1))
}

#' Determine feature enrichment per class using truth or predicted labels
#'
#' This function identifies the top N features (genes) for
#' each subtype based on their prevalence in the dataset using either
#' the truth labels or the predicted subgroups from DLBCLone.
#'
#' @param sample_metadata Data frame containing sample metadata with class labels,
#' by default in a column named "lymphgen". Use `label_column` to specify a different column.
#' @param label_column Name of the column containing the
#' class labels. The default is to use "lymphgen", the default truth class. 
#' @param truth_classes Vector of class labels to consider (default: c("BN2","EZB","MCD","ST2","N1")).
#' @param method Method to determine top features: "frequency" for most abundant
#' features, "chi_square" for top differentially mutated features in the classes
#' vs all other classes (default : "frequency").
#' @param num_feats Number of top features to display per subtype (default: 10).
#' @param separate_plot_per_group If TRUE, creates separate plots for each group
#' and also combines them with ggarrange (default: FALSE).
#' @param title Title for the plot (default: NULL).
#' @param p_threshold Maximum P value to retain (when method is fisher)
#' @param base_size Base font size used (passed to theme_Morons)
#' 
#'
#' @return A ggplot2 object representing the stacked bar plot.
#' 
#' @import dplyr ggplot2 tidyr rlang ggpubr
#' @export
#'
#' @examples
#' \dontrun{
#' library(GAMBLR.predict)
#' 
#' # Assuming my_DLBCLone_opt is the output from DLBCLone_optimize_params
#' plot_list <- posthoc_feature_enrichment(
#'     my_DLBCLone_opt$predictions,
#'     features=DLBCLone_model$features,
#'     method = "chi_square",
#'     num_feats = 10,
#'     title = "LymphGen"
#' ) 
#' print(plot_list$bar_plot)
#' 
#' plot_list <- posthoc_feature_enrichment(
#'    label_column = "lymphgen",
#'    sample_metadata = my_DLBCLone_opt$predictions,
#'    features = my_DLBCLone_opt$features,
#'    method = "fisher",
#'    num_feats = 10,
#'    base_size=9
#' )
#' print(plot_list$forest_plot)
#' }
#'
posthoc_feature_enrichment <- function(
  sample_metadata,
  features,
  label_column = "lymphgen",
  truth_classes = c("BN2","EZB","MCD","ST2","N1"),
  method = "frequency",
  num_feats = 10,
  p_threshold = 0.01,
  title = NULL,
  base_size = 7,
  separate_plot_per_group = TRUE
){


  bad_cols <- colSums(features) <= 0.02 * nrow(features)
  features <- features[, !bad_cols]
  features_binary = features %>% select(-ends_with("_feats")) #in case meta-features are present
  features_binary[features_binary>0]=1
  annotated_feats <- features_binary %>%
    rownames_to_column("sample_id") %>%
    left_join(
      sample_metadata %>% select(sample_id, !!sym(label_column)), 
      by = "sample_id"
    ) %>%
    filter(!is.na(!!sym(label_column)),!!sym(label_column) %in% truth_classes) %>%
    column_to_rownames("sample_id") 

  annotated_feats[[label_column]] <- as.factor(annotated_feats[[label_column]])
  
  gene_cols <- setdiff(colnames(annotated_feats), c("sample_id", label_column))
  subtypes <- truth_classes

  top_genes_per_subtype <- list()
  to_return = list()
  if(method == "frequency"){

    for(subtype in subtypes){
      subtype_samples <- annotated_feats[annotated_feats[[label_column]] == subtype, ]
      
      gene_counts <- colSums(subtype_samples[, gene_cols, drop = FALSE])
      gene_ranking <- order(gene_counts, decreasing = TRUE)
      
      # Top genes for this subtype
      top_genes_per_subtype[[subtype]] <- gene_cols[gene_ranking[1:num_feats]]
    }
  }else if(grepl("chi_square",method)){

    for(subtype in subtypes){
      # Create binary class: subtype vs rest
      y_binary <- factor(ifelse(annotated_feats[[label_column]] == subtype, subtype, paste0("not_", subtype)))
  
      chi_stats <- numeric(length(gene_cols))
  
      for(i in seq_along(gene_cols)){
        gene <- gene_cols[i]
        tbl <- table(annotated_feats[[gene]], y_binary)
        test <- suppressWarnings(chisq.test(tbl))
        chi_stats[i] <- test$`p.value`
      }
  
      gene_ranking <- order(chi_stats, decreasing = FALSE)
  
      # Top genes for this subtype
      top_genes_per_subtype[[subtype]] <- gene_cols[gene_ranking[1:num_feats]]
    }
    
  }else if(grepl("fisher",method)){
    inout_df = data.frame()
    annotated_feats_no_other = filter(annotated_feats,!!sym(label_column)!="Other")
    for(subtype in subtypes){
      # Create binary class: subtype vs rest
      y_binary <- factor(ifelse(annotated_feats_no_other[[label_column]] == subtype,
        subtype, paste0("not_", subtype)),levels=c(paste0("not_", subtype),subtype))
  
      for(i in seq_along(gene_cols)){
        gene <- gene_cols[i]
        tbl <- table(annotated_feats_no_other[[gene]], y_binary)
        
        test <- suppressWarnings(fisher.test(tbl))
        fisher.res = broom::tidy(test)
        fisher.res$group = subtype
        fisher.res$gene = gene
        inout_df = bind_rows(inout_df,fisher.res)
      }
    }
    p_threshold = 0.01
    inout_df = arrange(inout_df,desc(estimate)) %>% 
      filter(estimate > 1,p.value < p_threshold) %>%
      group_by(group) %>% 
      slice_head(n=num_feats) %>%
      ungroup()
    if(separate_plot_per_group){
        use_global_x <- TRUE        # Consider making this an optional argument. For now, leave hardcoded
        pad_mult     <- 0.01        # 1% padding on each side

        global_xlim <- NULL
        if (use_global_x) {
          df_all <- dplyr::filter(inout_df, group %in% subtypes)

          # collect all x values used across panels (on the *log* scale used in the plot)
          x_all <- c(log(df_all$estimate), log(df_all$conf.low), log(df_all$conf.high))
          x_all <- x_all[is.finite(x_all)]

          if (length(x_all)) {
            rng <- range(x_all)
            pad <- diff(rng) * pad_mult
            global_xlim <- c(rng[1] - pad, rng[2] + pad)
          }
        }
     # per-group plots
        plots <- list()

        cols <- get_gambl_colours()
        cols <- cols[match(truth_classes, names(cols))]
        names(cols) <- truth_classes

        last_idx <- length(subtypes)
        row_sizes <- numeric()
        max_n <- num_feats
        min_rows <- 3 
        overhead_rows <- 0
        cap_len <- grid::unit(0.5, "mm")

        for (i in seq_along(subtypes)) {
          subtype <- subtypes[[i]]

          sub_df <- inout_df |>
            dplyr::filter(group == subtype) |>
            dplyr::arrange(dplyr::desc(estimate)) |>
            dplyr::slice_head(n = max_n) |>
            dplyr::mutate(group = factor(group, levels = truth_classes))
          

          # order genes for this panel
          real_levels <- rev(unique(sub_df$gene))
          
          need_pad <- max(0, min_rows - length(real_levels))
          pad_levels <- if (need_pad > 0) paste0("<<pad_", seq_len(need_pad), ">>") else character(0)
          sub_df$gene <- factor(sub_df$gene, levels = c(real_levels, pad_levels))

          eff_rows <- max(length(real_levels), min_rows)
          row_sizes <- c(row_sizes, (eff_rows + overhead_rows) / (max_n + overhead_rows))

          # ----- legend seed: make sure x is always finite -----
          x_candidates <- c(log(sub_df$estimate), log(sub_df$conf.low), log(sub_df$conf.high))
          x_candidates <- x_candidates[is.finite(x_candidates)]
          x_dummy <- if (length(x_candidates)) mean(x_candidates) else 0  # safe fallback

          legend_seed <- data.frame(
            group = factor(truth_classes, levels = truth_classes),
            gene  = factor(rep(levels(sub_df$gene)[1], length(truth_classes)),
                          levels = levels(sub_df$gene)),
            x     = x_dummy
          )

          p_sub <- ggplot(sub_df, aes(x = log(estimate), y = gene, colour = group)) +
          geom_segment(
            aes(x = log(conf.low), xend = log(conf.high), y = gene, yend = gene),
            arrow = grid::arrow(length = cap_len, angle = 90, ends = "both", type = "open"),
            lineend = "butt",
            show.legend = FALSE   # <- so our segments don't contribute to legend graphic
          ) +
          geom_point() +
          # legend seeding layer you already have (keep as-is) ...
          geom_point(data = legend_seed, aes(x = x, y = gene, colour = group),
                    inherit.aes = FALSE, alpha = 0, show.legend = TRUE) +
          scale_colour_manual(
            name   = NULL,
            values = cols,
            limits = truth_classes,
            breaks = truth_classes,
            drop   = FALSE,
            na.translate = FALSE
          ) +
          guides(colour = guide_legend(
            title = NULL,
            ncol = 1,
            byrow = TRUE,                    # <- vertical legend
            override.aes = list(shape = 16, size = 3,  # <- point appearance in legend
                                linetype = 0, alpha = 1)
          )) +
            theme_Morons(base_size = base_size) +
            theme(legend.position = "none",axis.title.y = element_blank()) 

          p_sub <- p_sub +
          scale_y_discrete(
            limits = c(real_levels, pad_levels),
            breaks = real_levels,
            expand = ggplot2::expansion(add = c(0.5, 0.5))  # ~0.5 row below and above
          )
          if (i != last_idx) {
            p_sub <- p_sub +
              theme(
                axis.title.x = element_blank(),
                plot.margin  = ggplot2::margin(t = 2, r = 5, b = 0, l = 5, unit = "pt")
              )
          } else {
            p_sub <- p_sub +
              labs(x = "log(estimate)") +
              theme(plot.margin = ggplot2::margin(t = 2, r = 5, b = 2, l = 5, unit = "pt"))
          }
          if (!is.null(global_xlim)) {
            # same x range for every panel; prevents dropping data & keeps arrow caps visible
            p_sub <- p_sub +
              scale_x_continuous(expand = expansion(mult = 0)) +
              coord_cartesian(xlim = global_xlim, clip = "off")
          }
          plots[[i]] <- p_sub
        }

        p <- ggpubr::ggarrange(plotlist = plots, ncol = 1,
                              heights = row_sizes,
                              align = "v",
                              common.legend = TRUE,
                              legend = "right") + 
                              ggtitle(label= paste("Top", num_feats, "genes per subtype by", method))

    } else{

      p = ggplot(inout_df,aes(x=log(estimate),y=gene,colour=group)) + 
        geom_point() + 
        geom_segment(aes(x = log(estimate), y = gene,
                     xend = log(conf.high), yend = gene),
                     arrow = arrow(angle=90,length = unit(0.02, "npc"))) +
    geom_segment(aes(x = log(estimate), y = gene,
                     xend = log(conf.low), yend = gene),
                     arrow = arrow(angle=90,length = unit(0.02, "npc"))) +
        scale_colour_manual(values=get_gambl_colours()) +
        facet_wrap(~group,scales="free_y") + 
        theme_Morons(base_size = base_size)
    }


    to_return[["forest_plot"]] = p
    to_return[["results"]] = inout_df
    for(subtype in subtypes){
      top_genes_per_subtype[[subtype]] <- filter(inout_df,group==subtype) %>% pull(gene)
    }

  } else{
    stop(
      paste(
        "Must specify a valid method:",
        "'frequency'     - for most abundant features",
        "'fisher_test' - Use Fisher's exact test",
        "'chi_square' - for significance in-class vs out-of-class",
        sep = "\n"
      )
    )
  }
  
  plot_df <- bind_rows(lapply(names(top_genes_per_subtype), function(subtype){
    genes <- top_genes_per_subtype[[subtype]]
    test = annotated_feats[[label_column]] == subtype
    subtype_samples <- annotated_feats[annotated_feats[[label_column]] == subtype, ]
    n_subtype <- nrow(subtype_samples)  # total samples in this subtype
  
    counts <- colSums(subtype_samples[, genes, drop = FALSE] > 0)  # ensure 0/1
    plot_df <- tibble(
      subtype = subtype,
      gene = names(counts),
      count = as.numeric(counts),
      prop_gene = as.numeric(counts) / n_subtype
  )
  }))
  plot_df <- plot_df %>%
    group_by(subtype) %>%
    mutate(
      prop = count / sum(count),
      rank = rank(-count, ties.method = "first")  # rank genes within subtype
    ) %>%
    ungroup()

  # cumulative position for placing labels
  plot_df <- plot_df %>%
    group_by(subtype) %>%
    arrange(rank) %>%
    mutate(
      cum_prop = cumsum(prop),         
      pos = cum_prop - prop / 2       
    ) %>%
    ungroup()

  colours <- get_gambl_colours()
  n_colours <- max(plot_df$rank)

  fill_map <- plot_df %>%
    group_by(subtype) %>%
    mutate(
      base_col = ifelse(!is.na(colours[subtype]), colours[subtype], "grey"),
      ramp = list(colorRampPalette(c(first(base_col), "white"))(max(rank))),
      fill_color = ramp[[1]][rank]
    ) %>%
    ungroup()

  # merge colors back into plot_df
  plot_df$fill_color <- fill_map$fill_color

  p = ggplot(plot_df, aes(x = subtype, y = prop, fill = fill_color)) +
    geom_bar(stat = "identity", color = "black", position = position_stack(reverse = TRUE)) +
    geom_text(
      aes(label = paste0(gene, " ", scales::percent(prop_gene, accuracy = 1))),
      position = position_stack(vjust = 0.5, reverse = TRUE),
      size = 3
    ) +
    scale_fill_identity() +   
    labs(
      title = paste(title, " Top", num_feats, "genes per subtype (method:", method, ")"),
      x = "Subtype",
      y = "Proportion of Samples with Mutation"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 1),
      legend.position = "none"
    )
  to_return[["bar_plot"]] = p
  return(to_return)
}