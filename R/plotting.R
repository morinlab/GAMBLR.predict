#' Heatmap visualization of mutations in nearest neighbors for a sample
#'
#' Generates a heatmap of feature values for the nearest neighbors of a specified sample,
#' based on a DLBCLone model object. This visualization helps to inspect the feature profiles
#' of samples most similar to the query sample.
#'
#' @param this_sample_id Character. The sample ID for which to plot the nearest neighbor heatmap.
#' @param DLBCLone_model List. A DLBCLone model object, either from \code{DLBCLone_optimize_params}
#'   or \code{DLBCLone_KNN} (with \code{predict_unlabeled = TRUE}).
#'
#' @return A ComplexHeatmap object showing the feature matrix for the nearest neighbors of the sample.
#'
#' @details
#' - For models of type \code{DLBCLone_optimize_params}, uses the \code{neighbors} and \code{lyseq_status} fields.
#' - For models of type \code{DLBCLone_KNN}, uses the \code{unlabeled_neighbors} field.
#' - The function extracts the feature matrix rows corresponding to the nearest neighbors of the specified sample,
#'   and plots a heatmap of features with nonzero values.
#'
#' @importFrom dplyr filter
#' @import ComplexHeatmap
#' @importFrom grid gpar
#' @importFrom circlize colorRamp2
#'
#' @examples
#' # Assuming 'predicted_out' is a the output of DLBCLone_KNN_predict
#' and "SomeSample_ID" is a valid sample ID from among the test (not training) samples
#' \dontrun{
#' nearest_neighbor_heatmap("SomeSample_ID", predicted_out)
#'}
nearest_neighbor_heatmap <- function(this_sample_id,
                                     DLBCLone_model,
                                     truth_column = "lymphgen",
                                     clustering_distance = "binary",
                                     font_size = 14){
  pred_name = NULL
  if(!missing(DLBCLone_model) && "type" %in% names(DLBCLone_model)){
    if(DLBCLone_model$type == "DLBCLone_optimize_params"){
      neighbor_df = DLBCLone_model$neighbors
      lyseq_status = DLBCLone_model$lyseq_status
    }else if(DLBCLone_model$type == "DLBCLone_KNN"){
      if(!"unlabeled_neighbors" %in% names(DLBCLone_model)){
        print(names(DLBCLone_model))
        stop("DLBCLone_model must be the output of DLBCLone_KNN with predict_unlabeled = TRUE")
      }else if(is.null(DLBCLone_model$unlabeled_neighbors)){
        #no neighbors found for any of the incoming samples
        message("No neighbors found for any sample. Returning NULL.")
        return(NULL)
      }
      neighbor_df = DLBCLone_model$unlabeled_neighbors
      neighbor_transpose = filter(neighbor_df,sample_id==this_sample_id) 
      if(nrow(neighbor_transpose) == 0){
        message("No neighbors found for sample ", this_sample_id, ". Returning NULL.")
        return(NULL)
      }
      neighbor_transpose = neighbor_transpose %>% t()
      #deal with fewer than K neighbours
      neighbor_transpose = neighbor_transpose[!is.na(neighbor_transpose)]
      pred_name = "DLBCLone_ko"
    }

  }else{
    stop("DLBCLone_model must be the output of DLBCLone_optimize_params or DLBCLone_KNN")
  }
  
  xx=DLBCLone_model$features_df[neighbor_transpose,]
  if(any(is.na(rownames(xx)))){
    print(neighbor_transpose)
    stop("something went wrong. Some samples are missing from features_df")
  }
  top = max(xx)
  mid = top/2
  col_fun = circlize::colorRamp2(c(0, mid, top), c("white", "#FFB3B3", "red"))
  if(!truth_column %in% colnames(DLBCLone_model$predictions)){
    print(head(DLBCLone_model$predictions))
    stop("missing",truth_column)
  }
  row_df = select(DLBCLone_model$predictions, sample_id, !!sym(truth_column), !!sym(pred_name)) %>% 
    filter(sample_id %in% rownames(xx)) 
  if(!this_sample_id %in% row_df$sample_id){
    if("unlabeled_predictions" %in% names(DLBCLone_model)){
      row_df = bind_rows(row_df,
                        select(DLBCLone_model$unlabeled_predictions, sample_id, !!sym(truth_column), !!sym(pred_name)) %>% 
                            filter(sample_id %in% rownames(xx))
                        )
       sample_class = filter(DLBCLone_model$unlabeled_predictions, sample_id == this_sample_id) %>%
        pull(!!sym(pred_name))
             
    }else{
      sample_class = NULL
      row_df = bind_rows(row_df,
                       tibble(sample_id = this_sample_id,
                       !!rlang::sym(truth_column) := NA_character_,
                       !!rlang::sym(pred_name) := NA_character_)) 
    }
  }
  row_df = row_df %>% 
    column_to_rownames("sample_id") 
  anno_colours = get_gambl_colours()
  anno_list = list()

  anno_list[[truth_column]] = anno_colours
  if (!is.null(pred_name)) {
    anno_list[[pred_name]] = anno_colours
  }
  xx = xx[,colSums(xx)>0,drop=FALSE]
  row_anno = rowAnnotation(
    df = row_df[rownames(xx),,drop=FALSE],
    col = anno_list,
    #annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = font_size),
    show_legend = FALSE
  )

  title_text = paste("Sample", this_sample_id, "classified as", sample_class)
  Heatmap(xx,
          col = col_fun,
          right_annotation = row_anno,
          clustering_distance_rows = clustering_distance,
          show_heatmap_legend = FALSE,
          column_title = title_text,
          column_title_gp = gpar(fontsize=font_size),
          column_names_gp = gpar(fontsize=font_size))
}

#' Basic UMAP Scatterplot
#'
#' Generates a simple UMAP scatterplot for visualizing sample clustering or separation.
#'
#' @param optimized Data frame containing at least V1, V2, sample_id, and grouping columns.
#' @param plot_samples Optional character vector of sample_ids to label in the plot
#' @param colour_by Column name to color points by. Defaults to `truth_column`.
#' @param truth_column Name of the truth/ground-truth column (default: "lymphgen").
#' @param pred_column  Name of the predicted-class column (default: "DLBCLone_ko").
#' @param other_label  Label used for the outgroup/unclassified class (default: "Other").
#' @param title Plot title.
#' @param use_plotly Logical; if FALSE and `plot_samples` provided, draw static labels.
#' @param custom_colours Optional named vector of colors for groups; falls back to `get_gambl_colours()`.
#'
#' @return A ggplot object.
#' @export
#' 
#' @examples 
#' 
#' \dontrun{
#' my_umap = make_and_annotate_umap(my_data, my_metadata)
#' 
#' basic_umap_scatterplot(my_umap$df, #the data frame containing V1 and V2 from UMAP
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
#' @param optimized_model List. The output from DLBCLone optimization, containing predictions, features, and metadata.
#'
#' @details
#' - Creates a directory for results if it does not exist.
#' - Saves UMAP scatterplots for all samples and for non-"Other" samples.
#' - Generates alluvial plots for different DLBCLone predictions
#' - Exports an oncoplot summarizing mutation and classification results.
#' - Uses `make_umap_scatterplot`, `make_alluvial`, and `prettyOncoplot` for visualization.
#'
#' @return No return value. Side effect: writes multiple PDF files to disk.
#'
#' @import ggalluvial
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
#' Generates a UMAP plot highlighting the neighborhood of a given sample, showing its nearest neighbors and their group assignments.
#'
#' @param single_sample_prediction_output A list containing prediction results and annotation data frames. 
#'        Must include elements \code{prediction} (data frame with prediction results) and \code{anno_df} (data frame with UMAP coordinates and annotations).
#' @param training_predictions The equivalent data frame of prediction results for all training samples (e.g. optimized_model$df)
#' @param this_sample_id Character. The sample ID for which the neighborhood plot will be generated.
#' @param prediction_in_title Logical. If \code{TRUE}, includes the predicted label in the plot title.
#' @param add_circle Plot will include a circle surrounding the set of neighbors. Set to FALSE to disable.
#' @param label_column Does nothing, i.e. this is not currently working.
#'
#' @return A \code{ggplot2} object representing the UMAP plot with the selected sample and its neighbors highlighted.
#'
#' @details
#' The function extracts the nearest neighbors of the specified sample, draws segments connecting the sample to its neighbors, and colors points by group (e.g., lymphgen subtype). The plot title can optionally include the predicted label.
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom rlang sym
#'
#' @export
#' @examples
#' 
#' # Assuming 'optimization_result' is the output of DLBCLone_optimize_params
#' # and 'output' is the result of DLBCLone_predict_single_sample
#' # on sample_id "SAMPLE123":
#' \dontrun{
#'  make_neighborhood_plot(output, optimization_result$df, "SAMPLE123")
#' }
make_neighborhood_plot <- function(single_sample_prediction_output,
                                   training_predictions,
                                  this_sample_id,
                                  prediction_in_title = TRUE,
                                  add_circle = TRUE,
                                  label_column = "DLBCLone_io",
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
  if(missing(training_predictions)){
    #training_predictions = single_sample_prediction_output$anno_df
    training_predictions = 
    left_join(select(single_sample_prediction_output$anno_df,sample_id,lymphgen),
      select(single_sample_prediction_output$df,sample_id,V1,V2))
  }else if(missing(single_sample_prediction_output)){

    #Just plot the single sample in the context of the rest based on the optimization
    single_sample_prediction_output = list()
    single_sample_prediction_output[["prediction"]] = filter(training_predictions, sample_id==this_sample_id) 
     single_sample_prediction_output[["anno_df"]] = training_predictions
  }else{
    single_sample_prediction_output$prediction = filter(single_sample_prediction_output$prediction, sample_id==this_sample_id)
  
  }
xmin = min(training_predictions$V1, na.rm = TRUE)
  xmax = max(training_predictions$V1, na.rm = TRUE)
  ymin = min(training_predictions$V2, na.rm = TRUE)
  ymax = max(training_predictions$V2, na.rm = TRUE)
  #extract the sample_id for all the nearest neighbors with non-Other labels
  my_neighbours = filter(single_sample_prediction_output$prediction,
                         sample_id == this_sample_id) %>% 
                  pull(neighbor_id) %>% strsplit(.,",") %>% unlist()

  #set up links connecting each neighbor to the sample's point
  links_df = filter(training_predictions,sample_id %in% my_neighbours) %>% mutate(group=lymphgen)
  my_x = filter(single_sample_prediction_output$anno_out,
                sample_id==this_sample_id) %>% pull(V1)
  my_y = filter(single_sample_prediction_output$anno_out,
                sample_id==this_sample_id) %>% pull(V2)
  if(prediction_in_title){
    title = paste(this_sample_id,
                  pull(single_sample_prediction_output$prediction,
                       !!sym(label_column)))
    if(single_sample_prediction_output$prediction[[label_column]] == "Other" && single_sample_prediction_output$prediction$predicted_label !="Other"){
      title = paste(title,"(",single_sample_prediction_output$prediction$predicted_label,")")
    }

  }else{
    title = this_sample_id
  }
  links_df = mutate(links_df,my_x=my_x,my_y=my_y)
  links_df = links_df %>% select(V1,V2,my_x,my_y,group) %>% mutate(length = abs(V1-my_x)+abs(V2-my_y))
  
  
  pp=ggplot(mutate(training_predictions,group=lymphgen),
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
#' @param df 
#' @param drop_composite 
#' @param colour_by 
#' @param drop_other 
#' @param high_confidence 
#' @param custom_colours 
#' @param add_labels 
#'
#' @returns
#' 
#' @import ggside
#' @export
#'
#' @examples
make_umap_scatterplot = function(df,
                                 drop_composite = TRUE,
                                 colour_by="lymphgen",
                                 drop_other = FALSE,
                                 high_confidence = FALSE,
                                 custom_colours,
                                 add_labels = FALSE,
                                 title = NULL){
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
  if(add_labels){
    labels = group_by(df,!!sym(colour_by)) %>%
      summarise(median_x = median(V1),median_y = median(V2)) 
  }
  unique_lg = unique(df[[colour_by]])
  if(any(!unique_lg %in% names(cols))){
    missing = unique_lg[!unique_lg %in% names(cols)]
    print(paste("missing colour for:",paste(missing,collapse=",")))

  }
  p = ggplot(df,
             aes(x=V1,y=V2,colour=!!sym(colour_by),label=cohort)) + 
             geom_point(alpha=0) + 
             geom_point(data=df %>% filter(!!sym(colour_by) =="Other"),alpha=0.8) + 
             geom_point(data=df %>% filter(!!sym(colour_by) !="Other"),alpha=0.8) + 
    scale_colour_manual(values=cols) + 
    scale_fill_manual(values=cols) + 
    theme_Morons() + 
    guides(colour = guide_legend(nrow = 1)) + 
    xlim(xmin,xmax) + 
    ylim(ymin,ymax) 
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
#' Computes overall accuracy, balanced accuracy, and sensitivity for predicted vs. true class labels.
#' Optionally excludes samples assigned to the "Other" class from accuracy calculations.
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
#' result <- report_accuracy(predictions_df)
#' result$overall
#' result$per_class
#'
#' @export
report_accuracy <- function(predictions,
                            truth = "lymphgen",
                            pred = "DLBCLone_io",
                            per_group = FALSE,
                            metric = "accuracy",
                            verbose = FALSE) {
  all_classes <- unique(predictions[[truth]])
  no_other_pred <- filter(predictions, !!sym(truth) != "Other")
  if (verbose) {
    print(paste(nrow(no_other_pred), "non-Other samples were assigned a label"))
    print(paste(
      filter(no_other_pred, !!sym(pred) != "Other") %>% nrow(),
      "non-Other samples were assigned to a non-Other class"))
  }

  conf_matrix_no <- confusionMatrix(
    factor(no_other_pred[[pred]], levels = all_classes),
    factor(no_other_pred[[truth]], levels = all_classes)
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
  overall <- conf_matrix$overall[["Accuracy"]]

  return(list(
    no_other = acc_no,
    per_class = bal_acc,
    mean_balanced_accuracy = mean(bal_acc, na.rm = TRUE),
    per_class_sensitivity = sensitivity,
    overall = overall,
    confusion_matrix_no_other = conf_matrix_no,
    confusion_matrix = conf_matrix
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
#' @param add_accuracy_to_title Logical; if TRUE, adds accuracy/concordance rate to the plot title.
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
#' @param concordant_label_relative_pos Position for concordant labels: 0 (left), 0.5 (middle), or 1 (right).
#'
#' @return A ggplot2 object representing the alluvial plot.
#'
#' @details
#' - Visualizes flows between original and predicted classes, highlighting concordant and discordant assignments.
#' - Annotates concordance rate, per-group accuracy, and unclassified rate as specified.
#' - Supports flexible labeling, coloring, and axis ordering for publication-quality plots.
#'
#' @examples
#' # Example usage:
#' # make_alluvial(optimized_result)
#'
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
    hide_legend = TRUE
) {
  predictions = optimized$predictions
  if(is.null(group_order)){
    group_order = optimized$truth_classes
    #print("setting group order:")
  }
  

  if(!is.null(optimized$truth_column)){
    truth_column = optimized$truth_column
    truth_name = truth_column
  }
  if(!is.null(optimized$pred_column)){
    pred_column = optimized$pred_column
     pred_name = pred_column
  }
  
 

  if (accuracy_per_group) {
    accuracies <- report_accuracy(predictions, 
        truth = truth_column,
        pred = pred_column,
        per_group = accuracy_per_group)
  }
  if(count_excluded_as_other){
    excluded_meta = optimized$sample_metadata_no_features %>%
      mutate(!!pred_column := "Other")
    predictions = bind_rows(excluded_meta,predictions)
  }
  if("total_samples_available" %in% names(optimized)){
      full_denominator = optimized$total_samples_available
  }else{
    full_denominator = nrow(predictions)
  }

  xx <- predictions %>%
    rename(
      !!truth_name := !!sym(truth_column),
      !!pred_name := !!sym(pred_column)
    ) 
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
    title = paste0(title," Concordance (non-Other): ",pc_conc,"%")
  }
  bacc = round(accuracies$mean_balanced_accuracy,4)
  
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
      mutate(flow_label = paste0(left_group, " concordant: ", num," (",percent,"%",")"))
    
  }
  if(add_unclass_rate){
    right_other = filter(lodes_mean_right,right_group=="Other",axis==2) %>%
      summarize(total=sum(num)) %>%
      mutate(unclass= round(100 * total / full_denominator,1)) %>%
      mutate(label=paste0(total,"\n (",unclass,"%)")) %>%
               mutate(!!sym(truth_name) := "Other",!!sym(pred_name) := "Other")
    unclass = pull(right_other,unclass)
    title = paste0(title,", Classification rate: ",100-unclass,"%")
    title = paste0(title," Bal Acc: ",bacc)
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
      
    ) +

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
      
    ) +
    scale_x_discrete(limits = c(truth_name, pred_name), expand = c(.1, .05)) +
    scale_fill_manual(values = group_colours) +
    scale_colour_manual(values = group_colours,
                        aesthetics = c("color", "segment.color")) +
    labs(
      y = "Number of Samples",
      x = NULL
    ) +
    theme_minimal() +
    ggtitle(title)
  if(hide_legend){
    pp = pp +
      theme(legend.position = "none") 
   
  }
  if(add_unclass_rate){
    pp = pp + 
      geom_label_repel(data=left_other,
                       direction = "x",
                       x=1,y=20,
                 aes(fill=!!sym(truth_name),label=label)) +
      geom_label_repel(data=right_other,
                       direction = "x",
                       x=2,y=20,
                 aes(fill=!!sym(pred_name),label=label))
  }
  pp
}
