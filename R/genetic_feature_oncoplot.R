
#' UMAP Feature Plot
#'
#' Plot UMAP coordinates colored by subtype with an overlay highlighting samples
#' positive for selected features.
#'
#' @param dlbclone_umap_out List output from \code{make_and_annotate_umap()} containing
#'   \code{$df} with UMAP coordinates and \code{$features} with feature values.
#' @param features Character vector of feature names to plot. If NULL, uses the
#'   full set of features underlying the UMAP output.
#' @param alpha Numeric transparency for background points.
#' @param combine Logical; if TRUE, return a combined patchwork plot.
#' @param ncol Optional number of columns when \code{combine=TRUE}.
#'
#' @return A list of ggplot objects, or a combined patchwork plot when
#'   \code{combine=TRUE}.
#'
#' @export
umap_feature_plot <- function(dlbclone_umap_out,
                              features = NULL,
                              alpha = 0.5,
                              combine = FALSE,
                              ncol = NULL) {

  meta_no_features <- dlbclone_umap_out$df %>% 
    select(sample_id,all_of(c("lymphgen","V1","V2")))
  meta_with_features <- left_join(meta_no_features, 
    dlbclone_umap_out$features %>% 
    rownames_to_column(var = "sample_id"), by = "sample_id")
  lgcol <- get_gambl_colours()

  if (is.null(features)) {
    features <- intersect(
      colnames(dlbclone_umap_out$features),
      colnames(meta_with_features)
    )
  }

  plot_list <- vector("list", length(features))
  names(plot_list) <- features

  for (i in seq_along(features)) {
    feature <- features[i]

    df <- meta_with_features
    df$.feature_state <- ifelse(df[[feature]] > 0, "JS3", "no")

    p <- ggplot(df, aes(V1, V2)) +
      geom_point(aes(colour = lymphgen), alpha = alpha) +
      geom_point(
        data = df[df$.feature_state == "JS3", , drop = FALSE],
        aes(colour = .feature_state),
        size = 1.3,
        show.legend = FALSE
      ) +
      scale_colour_manual(values = lgcol) +
      theme_minimal() +
      ggtitle(feature) +
      theme(legend.position = "none")

    plot_list[[i]] <- p
  }

  if (combine) {
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      stop("combine=TRUE requires the 'patchwork' package. Install it or set combine=FALSE.")
    }

    # DimPlot-like behaviour: if ncol not provided, use a sensible default
    ncol_use <- ncol %||% ceiling(sqrt(length(plot_list)))

    return(patchwork::wrap_plots(plot_list, ncol = ncol_use))
  }

  return(plot_list)
}

#' Genetic Feature Oncoplot
#'
#' Render a ComplexHeatmap oncoplot from DLBCLone UMAP output or a feature matrix
#' with sample metadata.
#'
#' @param dlbclone_umap_out Optional list output from \code{make_and_annotate_umap()}.
#' @param features Optional feature matrix/data.frame (rows = samples, cols = features).
#' @param sort_by_annotation Column in metadata to sort samples by (default: \code{"V1"}).
#' @param top_annotations Optional data frame/matrix of per-sample annotations to show above the plot.
#' @param these_samples_metadata Optional metadata data frame (must include \code{sample_id}).
#' @param metadataColumns Optional vector of column names in metadata to include in the output. 
#' @param gene_order Optional vector specifying gene order (currently unused).
#' @param sample_order Optional vector specifying sample order.
#' @param gene_metadata Optional data frame describing genes (e.g., \code{feature}, \code{lymphgen}).
#' @param ... Passed to \code{ComplexHeatmap::Heatmap()}.
#'
#' @export
genetic_feature_oncoplot = function(
                        dlbclone_umap_out = NULL,
                        features,
                        sort_by_annotation = "lymphgen",
                        top_annotations = NULL,
                        these_samples_metadata = NULL,
                        metadataColumns = "lymphgen",
                        gene_order = NULL,
                        sample_order = NULL,
                        gene_metadata = NULL,
                        ...){
  # TODO: Make this function work with the objects stored in the DLBCLone UMAP outputs (i.e. make_and_annotate_umap)
  #  or a user-provided combination of a feature matrix and sample metadata
  gcol = get_gambl_colours()
  if(!is.null(dlbclone_umap_out)){
    #drop the meta-features if any are present (end with _feats)
    features = dlbclone_umap_out$features %>% 
      select(-ends_with("_feats"))
  
    #extract all the metadata from DLBCLone output (if provided)
    annotations = dlbclone_umap_out$df
  }else{
    if (is.null(these_samples_metadata) || is.null(features)) {
      stop("Must provide both these_samples_metadata and features if dlbclone_umap_out is NULL")
    }
    annotations = these_samples_metadata %>% 
      filter(sample_id %in% rownames(features))
  }

  if(!is.null(these_samples_metadata)){
    annotations = filter(annotations,
                         sample_id %in% these_samples_metadata$sample_id)
  }
  
  annotations = arrange(annotations,!!sym(sort_by_annotation))

  
  match_indices <- match(annotations$sample_id, rownames(features))

  # Reorder df1 using the generated indices
  features = features[match_indices, ]
  tfeatures = t(features)  
  if(!is.null(sample_order)){
    #check for different numbers of samples
    #sample_order = colnames(tfeatures)[match(sample_order,colnames(tfeatures))]
    sample_order_prune = intersect(sample_order,colnames(tfeatures))
    match_indices <- match(sample_order_prune, colnames(tfeatures))

    table(is.na(match_indices))
    # Reorder df1 using the generated indices
    tfeatures = tfeatures[,match_indices]
    
    if(!is.null(top_annotations)){
      top_annotations = top_annotations[sample_order_prune,]
    }
    
    if(!length(sample_order_prune)== ncol(tfeatures)){
      print(paste(length(sample_order),"!=",ncol(tfeatures)))
      stop()
    }
  }else{
    sample_order = arrange(annotations,!!sym(sort_by_annotation)) %>% pull(sample_id)
    sample_order_prune = intersect(sample_order,colnames(tfeatures))
    match_indices <- match(sample_order_prune, colnames(tfeatures))
    tfeatures = tfeatures[,match_indices]
    top_annotations = top_annotations[sample_order_prune,]
    #print(dim(tfeatures))
    #print(top_annotations)
  }

  if(!is.null(top_annotations)){
    top_ribbon <- HeatmapAnnotation(
      proportion =
        anno_barplot(
          top_annotations,
          gp = gpar(fill = gcol[colnames(top_annotations)],
                    col = NA),
          bar_width = 1,
          border = FALSE
        ),
      annotation_height = unit(12, "mm")
    )
  }else{
    top_ribbon = NULL
  }  
  ha_list = Heatmap(matrix(nrow = 0, ncol = ncol(tfeatures)), 
                    top_annotation = top_ribbon)
  #ha_list = NULL
  
  if(!is.null(gene_metadata)){
    unique_feats = unique(gene_metadata$lymphgen)
  }else{
    gene_metadata = data.frame(feature=rownames(tfeatures),lymphgen=NA)
    unique_feats = "unknown"
  }
  for(f in unique_feats){
    col_fun = colorRamp2(c(0, 1, 2), c("white", gcol[f],gcol[f]))
    
       f_rows = filter(gene_metadata, lymphgen== f) %>%
             pull(feature)
         
         lymphgen_heatmap = Heatmap(tfeatures[f_rows,,drop=FALSE],
                                    cluster_columns=FALSE,
                                    show_column_names=FALSE,
                                    row_names_gp =   gpar(fontsize=5),
                                    col = col_fun,
                                    ...
                                    )
    ha_list = ha_list %v%   lymphgen_heatmap
  }
  if(any(is.na(gene_metadata$lymphgen))){
    col_fun = colorRamp2(c(0, 1, 2), c("white", "grey", "orange"))
    
    f_rows = filter(gene_metadata, is.na(lymphgen)) %>%
      pull(feature)
    
    lymphgen_heatmap = Heatmap(tfeatures[f_rows,,drop=FALSE],
                               cluster_columns=FALSE,
                               show_column_names=FALSE,
                               row_names_gp =   gpar(fontsize=5),
                               col = col_fun,
                               ...
    )
    ha_list = ha_list %v%   lymphgen_heatmap
  }
  bottom_df = select(these_samples_metadata, sample_id, all_of(metadataColumns)) %>%
    column_to_rownames("sample_id")
  bottom_df = bottom_df[sample_order_prune,,drop=FALSE]
  bottom_anno = HeatmapAnnotation(
    df=bottom_df,
    col=list(lymphgen=gcol,DLBCLone_w=gcol,DLBCLone_wo=gcol)
    )
  ha_list = ha_list %v% Heatmap(matrix(nrow=0,
                                       ncol=ncol(tfeatures)),
                              bottom_annotation = bottom_anno)
  
  plot(ha_list,show_heatmap_legend=FALSE)

}
