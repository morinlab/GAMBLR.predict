
#' Summarize and Export DLBCLone Model Results
#'
#' Generates and saves a set of summary plots and tables for a DLBCLone model, including UMAP scatterplots, alluvial plots, and oncoplots.
#' Results are saved as PDF files in a directory named after the provided base name.
#'
#' @param base_name Character. The base name (and directory) for saving output files.
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
#' @examples
#' DLBCLone_summarize_model("Full_geneset_unweighted", optimized_model)
#'
#' @export
DLBCLone_summarize_model = function(base_name,optimized_model){
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
  mc = c("lymphgen","DLBCLone_wo","DLBCLone_io","DLBCLone_w","DLBCLone_i","DLBClass")
  sc = c("DLBCLone_wo","DLBCLone_io","lymphgen","DLBClass","DLBCLone_w","DLBCLone_i","confidence")
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
#' make_neighborhood_plot(output, optimization_result$df, "SAMPLE123")
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
  my_x = filter(single_sample_prediction_output$anno_df,
                sample_id==this_sample_id) %>% pull(V1)
  my_y = filter(single_sample_prediction_output$anno_df,
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
  if(!missing(custom_colours)){
    cols = custom_colours
  }else{
    cols = get_gambl_colours()
  }
  if(drop_composite){
    df = filter(df,!is.na(lymphgen),!grepl("COMP",lymphgen))
  }
  if(drop_other){
    df = filter(df,!is.na(lymphgen),lymphgen!="Other",lymphgen!="NOS")
  }
  if(high_confidence){
    df = filter(df,Confidence > 0.7)
  }
  if(add_labels){
    labels = group_by(df,!!sym(colour_by)) %>%
      summarise(median_x = median(V1),median_y = median(V2)) 
  }
  unique_lg = unique(df$lymphgen)
  if(any(!unique_lg %in% names(cols))){
    missing = unique_lg[!unique_lg %in% names(cols)]
    print(paste("missing colour for:",paste(missing,collapse=",")))

  }
  p = ggplot(df,
             aes(x=V1,y=V2,colour=!!sym(colour_by),label=cohort)) + geom_point(alpha=0.8) + 
    scale_colour_manual(values=cols) + theme_Morons() + 
    guides(colour = guide_legend(nrow = 1)) + 
    xlim(xmin,xmax) + 
    ylim(ymin,ymax) 
  if(!is.null(title)){
    p = p + ggtitle(title)
  }
  if(add_labels){
    p = p + geom_label_repel(data=labels,aes(x=median_x,y=median_y,label=!!sym(colour_by)))
  }
  ggMarginal(p,groupColour = TRUE,groupFill=TRUE)
  
}



#' Report Classification Accuracy and Per-Class Metrics
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
#'   \item{per_class}{Balanced accuracy per class}
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
report_accuracy = function(predictions,
                          truth="lymphgen",
                          pred="DLBCLone_io",
                          per_group = FALSE,
                          metric = "accuracy",
                          verbose = FALSE){
  if("BN2" %in% predictions[[pred]]){
    truth = "lymphgen"
  }
  all_classes = unique(predictions[[truth]])
  no_other_pred = filter(predictions,lymphgen != "Other")
  if(verbose){
    print(paste(nrow(no_other_pred),"non-Other samples were assigned a label"))
    print(paste(filter(no_other_pred,!!sym(pred) !="Other") %>% nrow(),"non-Other samples were assigned to a non-Other class"))
  }
  
  conf_matrix_no <- confusionMatrix(factor(no_other_pred[[pred]],levels=all_classes), 
    factor(no_other_pred[[truth]],levels=all_classes))
  #print(conf_matrix)
  overall = conf_matrix_no$overall[["Accuracy"]]

    if(metric == "accuracy"){
      #bal_acc_no <- conf_matrix$byClass[, "Balanced Accuracy"]
      acc_no = overall
    }else{
      stop("unsupported metric")
    }
  
  conf_matrix <- confusionMatrix(factor(predictions[[pred]],levels=unique(predictions[[truth]])), factor(predictions[[truth]],levels=unique(predictions[[truth]])))
    if(metric == "accuracy"){
      bal_acc <- conf_matrix$byClass[, "Balanced Accuracy"]
      sensitivity = conf_matrix$byClass[, "Sensitivity"]
    }else{
      stop("unsupported metric")
    }
  overall = conf_matrix$overall[["Accuracy"]]

  return(list(no_other=acc_no,
              per_class=bal_acc,
              mean_balanced_accuracy = mean(bal_acc,na.rm=TRUE),
              per_class_sensitivity = sensitivity,
              overall=overall,
              confusion_matrix_no_other=conf_matrix_no,
              confusion_matrix=conf_matrix))
}


#' Create an Alluvial Plot Comparing Original and Predicted Classifications
#'
#' This function generates a detailed alluvial plot to visualize the concordance and discordance between original (e.g., Lymphgen) and predicted (e.g., DLBCLone) class assignments for samples.
#' It supports annotation of concordance rates, per-group accuracy, unclassified rates, and flexible labeling and coloring options.
#'
#' @param optimized List containing prediction results and metadata, typically output from a DLBCLone optimization function.
#' @param pred Name of the column in predictions to use for the predicted class (default: "predicted_label_optimized").
#' @param count_excluded_as_other Logical; if TRUE, samples excluded due to missing features are counted as "Other" in the plot.
#' @param title Plot title (default: empty string).
#' @param group_order Character vector specifying the order of groups/classes for axes and coloring.
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
#' @export
make_alluvial <- function(
    optimized,
    count_excluded_as_other = FALSE,
    title = "",
    group_order = NULL,
    add_accuracy_to_title = TRUE,
    accuracy_per_group = TRUE,
    accuracy_type = "sensitivity",
    truth_name = "Lymphgen",
    truth_column = "lymphgen",
    pred_name = "DLBCLone",
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
    verbose = FALSE
) {
  predictions = optimized$predictions
  if(is.null(group_order)){
    group_order = optimized$truth_classes
    print("setting group order:")
  }
  
  print(group_order)

  if(is.null(truth_column) && !is.null(optimized$truth_column)){
    truth_column = optimized$truth_column
  }
  if (accuracy_per_group) {
    accuracies <- report_accuracy(predictions, 
        truth = truth_column,
        pred = pred_column,
        per_group = accuracy_per_group)
  }
  print("Done accuracyies")
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
    #group_by(lymphgen, predicted_label_optimized) %>%
    #summarize(num = n(), .groups = "drop") %>%
    rename(
      !!truth_name := !!sym(truth_column),
      #!!new_name := predicted_label_optimized
      !!pred_name := !!sym(pred_column)
    ) 
  print("Done renaming")
  xx <- xx %>%
    group_by(!!sym(truth_name), !!sym(pred_name)) %>%
    summarize(num = n(), .groups = "drop")


  grid <- expand_grid(
    !!sym(truth_name) := group_order,
    !!sym(pred_name) := group_order
  )
  print(colnames(grid))
  print(colnames(xx))
  xx <- grid %>%
    left_join(xx, by = c(truth_name, pred_name)) %>%
    mutate(num = replace_na(num, 0))
print("Done joining")
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
  title = paste0(title," Concordance (non-Other): ",pc_conc,"%")
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


  #print("HERE")
  lodes_mean_left = group_by(left_lodes,right_group,left_group) %>%
    mutate(y=mean(flow_y)) %>% ungroup() %>%
    arrange(right_group,left_group) %>%
    mutate(nudge_dir = -nudge)
    #print("HERE")
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
   #print("HERE")
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
#print(lodes_match)
  #print("HERE!")
xx_denom <- predictions %>%
  group_by(!!sym(truth_column)) %>%
  summarize(denom = n(), .groups = "drop") %>%
  rename(
    stratum = !!sym(truth_column)
  )
  #print("HERE!!") 
if(add_percent){
    lodes_match = left_join(lodes_match,xx_denom,by="stratum") %>%
      mutate(percent = round(100*num/denom,1)) %>%
      mutate(flow_label = paste0(flow_label," (",percent,"%",")"))
  }
  if(add_unclass_rate){
    right_other = filter(lodes_mean_right,right_group=="Other",axis==2) %>%
      summarize(total=sum(num)) %>%
      mutate(unclass= round(100 * total / full_denominator,1)) %>%
      mutate(label=paste0(total,"\n (",unclass,"%)")) %>%
               mutate(!!sym(truth_name) := "Other",!!sym(pred_name) := "Other")
    #print(right_other)
    unclass = pull(right_other,unclass)
    title = paste0(title,", Classification rate: ",100-unclass,"%")
    title = paste0(title," Bal Acc: ",bacc)
    left_other = filter(lodes_mean_left,left_group=="Other",axis==2) %>%
      summarize(total=sum(num)) %>%
      mutate(unclass= round(100 * total / full_denominator,1)) %>%
      mutate(label=paste0(total,"\n (",unclass,"%)")) %>%
      mutate(!!sym(truth_name) := "Other",!!sym(pred_name) := "Other")
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
    scale_fill_manual(values = get_gambl_colours()) +
    scale_colour_manual(values = get_gambl_colours(),
                        aesthetics = c("color", "segment.color")) +
    labs(
      y = "Number of Samples",
      x = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "none") + 
    ggtitle(title)
  
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