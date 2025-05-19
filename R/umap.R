#' Optimize the threshold for classifying samples as "Other"
#'
#' Performs a post-hoc evaluation of the classification of a sample as one of
#' the main classes vs the outgroup/unclassified label "Other" and returns the
#' optimal threshold for classifying a sample as "Other" based on the ground
#' truth provided in the true_labels vector. It evaluates the performance
#' of the classifier using a range of thresholds and returns the best threshold
#' based on the specified metric (balanced accuracy or accuracy). 
#' This function is not generally meant to be called directly but rather is
#' a helper function used by DLBCLone_optimize_params.
#'
#' @param predicted_labels Vector of predicted labels for the samples
#' @param true_labels Vector of true labels for the samples
#' @param other_score Vector of scores for the "Other" class for each sample ()
#' @param all_classes Vector of classes to use for training and testing.
#' Default: c("MCD","EZB","BN2","N1","ST2","Other")
#' @param maximize Metric to use for optimization.
#' Either "accuracy" (actual accuracy across all samples) or "balanced_accuracy"
#' (the mean of the balanced accuracy values across all classes).
#' Default: "balanced_accuracy"
#' @param exclude_other_for_accuracy Set to TRUE to exclude the
#' "Other" class from the 'lymphgen' column when calculating accuracy metrics
#' (passed to DLBCLone_optimize_params). Default: FALSE
#'
#' @returns a list of data frames with the predictions and the UMAP input
#' @export
#'
optimize_outgroup <- function(predicted_labels,
                            true_labels,
                            other_score,
                            all_classes = c("MCD",
                                            "EZB",
                                            "BN2",
                                            "N1",
                                            "ST2",
                                            "Other"),
                            maximize ="balanced_accuracy",
                            exclude_other_for_accuracy = FALSE){
  
  rel_thresholds = seq(1,10,0.1)
  sens_df = data.frame()
  acc_df = data.frame()
  predictions = data.frame(predicted_label=as.character(predicted_labels),
                           true_label=as.character(true_labels))

  for(threshold in rel_thresholds){
      predictions_new = mutate(predictions,
                               predicted_label = ifelse(other_score < threshold,
                                                        predicted_label,
                                                        "Other"))

      pred = factor(predictions_new[["predicted_label"]],levels=all_classes)
      truth = factor(predictions_new[["true_label"]],levels=all_classes)
      conf_matrix <- confusionMatrix(pred, truth)

      bal_acc <- conf_matrix$byClass[, "Balanced Accuracy"]
      if(maximize == "balanced_accuracy"){
        bal_acc$average_accuracy = mean(bal_acc)
      }else{
        bal_acc$average_accuracy = conf_matrix$overall[["Accuracy"]]
      }
      bal_acc$threshold = threshold
      acc_df = bind_rows(acc_df,bal_acc)
      sn <- conf_matrix$byClass[, "Sensitivity"]  
      sn$average_sensitivity = mean(sn)
      sn$threshold = threshold
      sens_df = bind_rows(sens_df,sn)
  }
  if(maximize %in% c("balanced_accuracy","accuracy")){
    best = slice_head(arrange(acc_df,desc(average_accuracy)),n=1)
  }else{
    best = slice_head(arrange(sens_df,desc(average_sensitivity)),n=1)

  }
  
  return(best)
}


#' Plot the result of a DLBCLone classification
#'
#' @param test_df Data frame containing the test data with UMAP coordinates
#' @param train_df Data frame containing the training data with UMAP coordinates
#' @param predictions_df Data frame containing the predictions with UMAP coordinates
#' @param other_df Data frame containing the predictions for samples in the "Other" class
#' @param details Single-row data frame with the best parameters from DLBCLone_optimize_params
#' @param annotate_accuracy Set to true to add labels with accuracy values
#' @param classes Vector of classes that were used in the training and testing
#' @param label_offset Length of the label offset for the accuracy labels
#' @param title1 additional argument
#' @param title2 additional argument
#' @param title3 additional argument
#'
#' @returns a ggplot object
#' @export
#'
#' @examples
#' #add the dataset name to the metadata if it's not already there (required for the plot work)
#' lymphgen_A53_DLBCLone$df$dataset = "GAMBL"
#'
#' DLBCLone_train_test_plot(
#'  test_df = lymphgen_A53_DLBCLone$df,
#'  train_df = lymphgen_A53_DLBCLone$df,
#'  predictions_df = lymphgen_A53_DLBCLone$predictions,
#'  #other_df = lymphgen_A53_DLBCLone$predictions_other, #required only when "Other" was in the truth_classes
#'  details = lymphgen_A53_DLBCLone$best_params,
#'  classes = c("MCD","EZB","BN2","ST2","N1","A53","Other"),
#'  annotate_accuracy=TRUE,label_offset = 1)
#'
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
      #print(details)
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


#' Run UMAP and attach result to metadata
#'
#' @param df Feature matrix with one row per sample and one column per mutation
#' @param metadata Metadata data frame with one row per sample and a column sample_id that
#' matches the row names of df. This data frame will be joined to the UMAP output.
#' @param umap_out Optional UMAP output from a previous run. If provided, the function
#' will use this model to project the data instead of re-running UMAP. This is useful
#' for reproducibility and for using the same UMAP model on different datasets.
#' @param join_column The column name in the metadata data frame that contains the sample IDs (default sample_id).
#' @param n_neighbors Passed to UMAP2. The number of neighbors to consider when calculating the UMAP embedding.
#' @param min_dist Passed to UMAP2. The minimum distance between points in the UMAP embedding.
#' @param metric Passed to UMAP2. The distance metric to use for calculating distances between points.
#' @param n_epochs Passed to UMAP2. The number of epochs to run the UMAP algorithm.
#' @param init Passed to UMAP2. The initialization method for the UMAP algorithm.
#' @param na_vals How to deal with NA values. Two options are "drop", which
#' will remove all columns containing at least one NA or "to_zero", which sets
#' all NA to zero and leaves the column intact.
#' @param seed Passed to UMAP2. The random seed for reproducibility.
#' @param ret_model additional argument
#'
#' @import uwot
#' @export
#'
#' @examples
#'
#' umap_outs = make_and_annotate_umap(df=gambl_train_lymphgen,
#'                            min_dist = 0,
#'                            n_neighbors = 55,
#'                            init="spectral",
#'                            n_epochs = 1500,
#'                            #seed=best_params$seed,
#'                            metadata=gambl_train_meta_dlbclass,
#'                            ret_model=T,
#'                            metric="cosine")
#'
#'
make_and_annotate_umap = function(df,
                              metadata,
                              umap_out,
                              n_neighbors=55,
                              min_dist=0,
                              metric="cosine",
                              n_epochs=1500,
                              init="spca",
                              ret_model=TRUE,
                              na_vals = "drop",
                              join_column="sample_id",
                              seed=12345,
                              target_column,
                              target_metric="euclidean",
                              target_weight=0.5){
  original_n = nrow(df)
  if(na_vals == "to_zero"){
    df[is.na(df)] = 0
  }else if(na_vals == "drop"){
    df <- df[, colSums(is.na(df)) == 0]
  }
  rs = rowSums(df,na.rm=TRUE)
  df = df[rs>0,]
  if(missing(df)){
    stop("provide a data frame or matrix with one row for each sample and a numeric column for each mutation feature")
  }
  if(missing(metadata)){
    stop("metadata is required and should contain a column sample_id that matches the row names of your mutation data frame")
  }
  keep_rows = rownames(df)[rownames(df) %in% metadata[[join_column]]]
  df= df[keep_rows,]
  metadata= filter(metadata,!!sym(join_column) %in% rownames(df))
  message(paste("kept",nrow(metadata),"rows of the data"))
  if(missing(umap_out)){
    if(missing(target_column)){
      umap_out = umap2(df %>% as.matrix(),
                       n_neighbors = n_neighbors,
                       min_dist = min_dist,
                       metric = metric,
                       ret_model = ret_model,
                       n_epochs=n_epochs,
                       init=init,
                       seed = seed,
                       n_threads = 1) # possibly add rng_type = "deterministic"
      #IMPORTANT: n_threads must not be changed because it will break reproducibility  
    }else{
      #supervised
      metadata[[target_column]] = factor(metadata[[target_column]])
      print(table(metadata[[target_column]]))
      umap_out = umap2(df %>% as.matrix(),
                       n_neighbors = n_neighbors,
                       min_dist = min_dist,
                       metric = metric,
                       ret_model = ret_model,
                       n_epochs=n_epochs,
                       init=init,
                       seed = seed,
                       n_threads = 1,
                       y = metadata[[target_column]],
                       target_metric = target_metric,
                       target_weight = target_weight
                       ) # possibly add rng_type = "deterministic"
      #IMPORTANT: n_threads must not be changed because it will break reproducibility
    }
    

  }else{
    umap_out = umap_transform(X=df,
                                    model=umap_out$model)
    ret_model = FALSE
  }
  if(ret_model){
    umap_df = as.data.frame(umap_out$embedding) %>% rownames_to_column(join_column)
  }else{
    umap_df = as.data.frame(umap_out) %>% rownames_to_column(join_column)
  }
  umap_df = left_join(umap_df,metadata)

  results = list()
  results[["df"]]=umap_df
  results[["features"]] = df
  if(ret_model){
    results[["model"]]= umap_out
  }
  return(results)
}

#' Optimize parameters for classifying samples using UMAP and k-nearest neighbor
#'
#' @param combined_mutation_status_df Data frame with one row per sample and
#' one column per mutation
#' @param metadata_df Data frame of metadata with one row per sample and
#' three required columns: sample_id, dataset and lymphgen
#' @param truth_classes Vector of classes to use for training and testing.
#' Default: c("EZB","MCD","ST2","N1","BN2","Other")
#' @param eval_group Specify whether certain rows will be evaluated and
#' held out from training rather than using all samples.
#' @param umap_out The output of a previous run of make_and_annotate_umap.
#' If provided, the function will use this model to project the data
#' instead of re-running UMAP.
#' @param min_k Starting k for knn (Default: 3)
#' @param max_k Ending k for knn (Default: 33)
#' @param verbose Whether to print verbose outputs to console
#' @param seed Random seed to use for reproducibility (default: 12345)
#' @param maximize Metric to use for optimization. Either "sensitivity"
#' (average sensitivity across all classes), "accuracy"
#' (actual accuracy across all samples) or "balanced_accuracy" (the mean of the
#' balanced accuracy values across all classes). Default: "balanced_accuracy"
#'
#' @returns List of data frames with the results of the parameter optimization
#' including the best model, the associated knn parameters and the annotated
#' UMAP output as a data frame. The list also includes the predictions for the
#' "Other" class if it was included in the training and testing.
#' 
#' @export
#'
#' @examples
#'
#' lymphgen_A53_DLBCLone =  DLBCLone_optimize_params(
#'    lgen_feat_status, #our binary feature matrix
#'    a53_meta, #our metadata
#'    umap_out = lymphgen_A53_all_feat_gambl, # force use existing UMAP fit
#'    eval_group = NULL, # use all samples for evaluating accuracy
#'    truth_classes = c("MCD","EZB","BN2","ST2","N1","A53","Other"))

DLBCLone_optimize_params = function(combined_mutation_status_df,
                           metadata_df,
                           umap_out,
                           truth_classes = c("EZB",
                                             "MCD",
                                             "ST2",
                                             "N1",
                                             "BN2",
                                             "Other"),
                           optimize_for_other = TRUE,
                           eval_group = "Lacy",
                           min_k=3,
                           max_k=33,
                           verbose = FALSE,
                           seed = 12345,
                           maximize = "balanced_accuracy",
                           exclude_other_for_accuracy = FALSE
                           ) {
  if(optimize_for_other){
    exclude_other_for_accuracy = FALSE
  }else{
    exclude_other_for_accuracy = TRUE
  }
  na_opt = c("drop")
  num_class = length(truth_classes)
  weights_opt = c(TRUE,FALSE)
  threshs = seq(0,0.9,0.1)
  ks = seq(min_k,max_k,2)
  results <- data.frame()
  best_params <- data.frame()

  use_w = TRUE
  best_acc = 0
  best_fit = NULL
  this_accuracy = 0
  best_pred = NULL
  other_pred = NULL
  for(na_option in na_opt){
    if(missing(umap_out)){
      outs = make_and_annotate_umap(df=combined_mutation_status_df,
                              min_dist = 0,
                              n_neighbors = 55,
                              n_epochs = 1500,
                              seed=seed,
                              metadata=metadata_df,
                              ret_model=T,
                              metric="cosine",
                              join_column="sample_id",
                              na_vals = na_option)
    }else{
      #project onto existing model instead of re-running UMAP
      outs = make_and_annotate_umap(df=combined_mutation_status_df,
                              umap_out = umap_out,
                              min_dist = 0,
                              n_neighbors = 55,
                              n_epochs = 1500,
                              seed=seed,
                              metadata=metadata_df,
                              ret_model=T,
                              metric="cosine",
                              join_column="sample_id",
                              na_vals = na_option)
    }
    ignore_top = FALSE
    if(is.null(eval_group)){
      ignore_top = TRUE
    }
    for(use_w in weights_opt){
      for(k in ks){
        message(paste("K:",k))
        for (threshold in threshs){
          if(is.null(eval_group)){
            test_coords = filter(outs$df,lymphgen %in% truth_classes) %>% select(V1,V2)
            train_coords = filter(outs$df,lymphgen %in% truth_classes) %>% select(V1,V2)
            train_labels = filter(outs$df,lymphgen %in% truth_classes) %>% pull(lymphgen)

          }else{

            test_coords = filter(outs$df,dataset == eval_group) %>% select(V1,V2)
            train_coords = filter(outs$df,dataset != eval_group,lymphgen %in% truth_classes) %>% select(V1,V2)
            train_labels = filter(outs$df,dataset != eval_group,lymphgen %in% truth_classes) %>% pull(lymphgen)
          }
          if(verbose){
            print(paste("k:",k,"threshold:",threshold,"use_weights:",use_w,"na_option:",na_option))
          }
          
          pred = weighted_knn_predict_with_conf(
            train_coords = train_coords,
            train_labels = train_labels,
            test_coords = test_coords,
            k=k,
            conf_threshold =threshold,
            na_label="Other",
            use_weights = use_w,
            ignore_top = ignore_top,
            verbose = verbose)

          if(is.null(eval_group)){
            xx_d = bind_cols(filter(outs$df,lymphgen %in% truth_classes) ,pred)
            train_d = filter(outs$df,lymphgen %in% truth_classes)
          }else{
            xx_d = bind_cols(filter(outs$df,dataset == eval_group) ,pred)
            train_d = filter(outs$df,dataset != eval_group,lymphgen %in% truth_classes)
          }
          if("Other" %in% truth_classes){
            xx_d$lymphgen = factor(xx_d$lymphgen)
          }else{
            if(is.null(eval_group)){
              test_coords = filter(outs$df,lymphgen %in% "Other") %>% select(V1,V2)
              train_coords = filter(outs$df,lymphgen %in% truth_classes) %>% select(V1,V2)
              train_labels = filter(outs$df,lymphgen %in% truth_classes) %>% pull(lymphgen)
            }else{
              test_coords = filter(outs$df,dataset == eval_group,lymphgen %in% "Other") %>% select(V1,V2)
              train_coords = filter(outs$df,dataset == eval_group,lymphgen %in% unique(c("Other",truth_classes))) %>% select(V1,V2)
              train_labels = filter(outs$df,dataset == eval_group,lymphgen %in% unique(c("Other",truth_classes))) %>% pull(lymphgen)
            }
            n_other = nrow(test_coords)

            if(!"Other" %in% truth_classes && n_other > 0){
              pred_other = weighted_knn_predict_with_conf(
                train_coords = train_coords,
                train_labels = train_labels,
                test_coords = test_coords,
                k=k,
                conf_threshold =threshold,
                na_label="Other",
                use_weights = use_w,
                ignore_top = ignore_top,
                verbose = verbose)

              if(!is.null(eval_group)){
                xx_o = bind_cols(filter(outs$df,dataset == eval_group,lymphgen =="Other") ,pred_other)
              }else{
                xx_o = bind_cols(filter(outs$df,lymphgen == "Other") ,pred_other)
              }
            }
            xx_d$lymphgen = factor(xx_d$lymphgen,levels = c(unique(xx_d$lymphgen),"Other"))

          }
          true_factor = xx_d$lymphgen
          pred_factor = factor(xx_d$predicted_label,levels = levels(xx_d$lymphgen))
          
          if(verbose){
            print("true_factor")
            print(levels(true_factor ))
            print("===")
            print(levels(xx_d$predicted_label))
          }

          conf_matrix <- confusionMatrix(pred_factor, true_factor)

          bal_acc <- conf_matrix$byClass[, "Balanced Accuracy"]  # one per class
          sn <- conf_matrix$byClass[, "Sensitivity"]  # one per class
          if(verbose){
            print(bal_acc)
          }

          overall_accuracy <- conf_matrix$overall[["Accuracy"]]
          if(exclude_other_for_accuracy){
            mean_balanced_accuracy = mean(bal_acc[!names(bal_acc) == "Class: Other"], na.rm = TRUE)
          }else{
            mean_balanced_accuracy = mean(bal_acc)
          }
          
          overall_sensitivity<- mean(sn[!names(sn) == "Class: Other"], na.rm = TRUE)
          
          if(optimize_for_other){
            optimized_accuracy_and_thresh = optimize_outgroup(pred_factor,
                                             true_factor,
                                             xx_d$other_score,
                                             all_classes = truth_classes,
                                             maximize = maximize,
                                             exclude_other_for_accuracy = exclude_other_for_accuracy)
            out_opt_thresh = optimized_accuracy_and_thresh$threshold
            optimized_accuracy_and_thresh$average_accuracy[is.na(optimized_accuracy_and_thresh$average_accuracy)] = 0
            out_opt_acc = optimized_accuracy_and_thresh$average_accuracy
          }else{
            out_opt_acc = 0
            out_opt_thresh = 0
          }
          if(maximize == "sensitivity"){
            this_accuracy = overall_sensitivity
          }else if(maximize == "balanced_accuracy"){
            this_accuracy =mean_balanced_accuracy
          }else{
            this_accuracy = overall_accuracy
          }
          if(out_opt_acc > this_accuracy){
            this_accuracy = out_opt_acc
          }
          row <- data.frame(k = k,
                            threshold = threshold,
                            use_weights = use_w,
                            optimized_accuracy = this_accuracy,
                            overall_accuracy = overall_accuracy,
                            mean_balanced_accuracy = mean_balanced_accuracy,
                            sensitivity = overall_sensitivity,
                            N1_sn = unname(sn["Class: N1"]),
                            BN2_sn= unname(sn["Class: BN2"]),
                            ST2_sn = unname(sn["Class: ST2"]),
                            BN2_bacc = unname(bal_acc["Class: BN2"]),
                            MCD_bacc = unname(bal_acc["Class: MCD"]),
                            EZB_bacc = unname(bal_acc["Class: EZB"]),
                            A53_bacc = unname(bal_acc["Class: A53"]),
                            Other_bacc = unname(bal_acc["Class: Other"]),
                            C1_bacc = unname(bal_acc["Class: C1"]),
                            C2_bacc = unname(bal_acc["Class: C2"]),
                            C3_bacc = unname(bal_acc["Class: C3"]),
                            C4_bacc = unname(bal_acc["Class: C4"]),
                            C5_bacc = unname(bal_acc["Class: C5"]),
                            ST2_bacc = unname(bal_acc["Class: ST2"]),
                            threshold_outgroup = out_opt_thresh, 
                            accuracy_out = out_opt_acc,
                            na_option= na_option) 
   
          
          if(this_accuracy > best_acc){
            best_acc = this_accuracy
   
            message(paste("best accuracy:",best_acc,"k:",k,"threshold:",threshold,"Other_thresh:",out_opt_thresh,"na:",na_option,"Balanced accuracy:",overall_accuracy, "sensitivity:",overall_sensitivity))

            best_fit = outs
            best_pred = xx_d
            if(!"Other" %in% truth_classes  && n_other > 0){
              other_pred = pred_other
            }

            best_params = row
          }
          results <- rbind(results, row)
        }
      }
    }
  }
  best_params$num_classes = num_class
  best_params$num_features = ncol(best_fit$features)
  best_params$seed = seed

  test_coords = outs$df %>% select(V1,V2)
  train_coords = filter(outs$df,lymphgen %in% truth_classes) %>% select(V1,V2)
  train_labels = filter(outs$df,lymphgen %in% truth_classes) %>% pull(lymphgen)
  pred = weighted_knn_predict_with_conf(
            train_coords = train_coords,
            train_labels = train_labels,
            test_coords = test_coords,
            k=best_params$k,
            conf_threshold =best_params$threshold,
            na_label="Other",
            use_weights = best_params$use_weights,
            ignore_top = ignore_top,
            verbose = verbose,
            track_neighbors = TRUE)
  xx_d = bind_cols(outs$df,pred)
  to_ret = list(params=results,
                best_params = best_params,
                model=best_fit$model,
                features=best_fit$features,
                df=outs$df,
                predictions=xx_d)
  if(!"Other" %in% truth_classes && n_other > 0){
    to_ret[["predictions_other"]] = xx_o
    to_ret[["predictions_combined"]] = bind_rows(xx_o,best_pred)
  }
  return(to_ret)
}



#' Weighted k-nearest neighbor with confidence estimate
#'
#' @param train_coords Data frame containing coordinates for samples with known
#' labels (training samples). One row per sample, one column per feature. Usually
#' this would be the two UMAP components
#' @param train_labels Vector of labels for training cases
#' @param test_coords Data frame containing coordinates for samples to be classified
#' with the same format as train_coords
#' @param k The number of neigbors to consider
#' @param epsilon This value is added to the distances before applying weights
#' when use_weights is TRUE. Default: 1e-5
#' @param conf_threshold Minimum confidence for classifying a sample based on
#' it's neighbors. Below this value the sample will be assigned na_label instead
#' @param na_label Class to assign all samples that are not confidently assigned.
#' Default: Other
#' @param verbose Whether to print verbose outputs to console
#' @param use_weights Set to FALSE for all neigbors to have equal weight when
#' calculating the confidence
#' @param ignore_top Set to TRUE to avoid considering a nearest neighbor with
#' distance = 0. This is usually only relevant when re-classifying labeled
#' samples to estimate overall accuracy
#' @param track_neighbors Set to TRUE to include details for the nearest neighbors of each sample
#'
#' @returns data frame with labels and confidence values for rows in test_coords
#' @export
#'
weighted_knn_predict_with_conf <- function(train_coords,
                                           train_labels,
                                           test_coords,
                                           k,
                                           epsilon = 0.1,
                                           conf_threshold = NULL,
                                           na_label = "Other",
                                           verbose = FALSE,
                                           use_weights = TRUE,
                                           ignore_top = FALSE,
                                           track_neighbors = TRUE,
                                           separate_other = TRUE) { #big change here. Other is considered separately for optimization
  if (nrow(train_coords)==0 || nrow(test_coords) == 0) {
    print("train_coords:")
    print(nrow(train_coords))
    print("test:")
    print(nrow(test_coords))
    stop("train_coords and test_coords must be data frames with at least one row")
  }
  # get the 100 nearest neighbors
  nn <- get.knnx(train_coords, test_coords, 100)
  all_neighbors = data.frame()
  preds <- character(nrow(test_coords))
  confs <- numeric(nrow(test_coords))

  train_labels = as.character(train_labels)
  for (i in 1:nrow(test_coords)) {
    if(verbose){
      print(paste("index:",i))
      print(test_coords[i,])
    }
    neighbors <- nn$nn.index[i, ]
    distances <- nn$nn.dist[i, ]
    if(ignore_top){
      #ignore a neighbor if it has identical V1 and V2
      distances <- nn$nn.dist[i, ]
      if(distances[1] == 0){
        neighbors = neighbors[-1]
        distances = distances[-1]
      }

    }
    distances = distances +  epsilon
    weights <- 1 / distances
    if(use_weights){
      weights <- 1 / distances
    }else{
      weights <- rep(1,length(distances))
    }


    neighbor_labels <- train_labels[neighbors]

    if(verbose){
      print("neighbors:")
      print(neighbors)
      print("distances:")
      print(distances)
      print("weights:")
      print(weights)
      print("labels:")
      print(neighbor_labels)
    }
    # Remove NAs (just in case)
    valid <- !is.na(neighbor_labels)
    #num_other_neighbors = sum(neighbor_labels == "Other")
    other_dists = distances[neighbor_labels == "Other"]
    if(separate_other){
      valid = valid & neighbor_labels != "Other"
    }
    neighbor_labels <- neighbor_labels[valid]
    weights <- weights[valid]
    distances <- distances[valid]
    neighbors <- neighbors[valid]
    #now take the first k neighbors
    if(length(neighbor_labels) > k){
      neighbor_labels = neighbor_labels[1:k]
      weights = weights[1:k]
      distances = distances[1:k]
      neighbors = neighbors[1:k]
    }

    others_closer = which(other_dists < max(distances))
    
    others_distances = other_dists[others_closer]
    if(use_weights){
      others_weights = 1 / others_distances
    }else{
      others_weights = rep(1,length(others_closer))
    }
    neighbors_other = length(others_closer)
    other_weighted_votes = sum(others_weights)
    #other_weighted_votes = neighbors_other
    #  print(paste("other weighted votes:",other_weighted_votes))
    #}
    if (length(neighbor_labels) == 0) {
      preds[i] <- "Other"
      confs[i] <- 1
      if(track_neighbors){
        
        rel_other = 10
        neighbor_info <- data.frame(
          other_score = rel_other,
          neighbor = paste(neighbors,collapse=","),
          distance = paste(round(distances, 3),collapse=","),
          label = paste(neighbor_labels,collapse=","),
          weighted_votes = "",
          neighbors_other = neighbors_other,
          other_weighted_votes = 0,
          total_w = 1,
          pred_w = 2
          
        )
        all_neighbors = bind_rows(all_neighbors, neighbor_info)
        next;
      }
    }

    weighted_votes <- tapply(weights, neighbor_labels, sum)
    
    if (length(weighted_votes) == 0) {
      preds[i] <- na_label
      confs[i] <- NA
      total_weight  <- 0
      pred_weight <- 0
    } else {
      predicted_label <- names(which.max(weighted_votes))
      total_weight <- sum(weighted_votes)
      pred_weight <- weighted_votes[predicted_label]
      confidence <- pred_weight / total_weight

      # Confidence thresholding
      if (!is.null(conf_threshold) && confidence < conf_threshold) {
        preds[i] <- na_label
        confs[i] <- confidence
      } else {
        preds[i] <- predicted_label
        confs[i] <- confidence
      }
    }

    if(track_neighbors){
      # Create a data frame to store neighbors, distances, and weights
      if(separate_other){
        rel_other = other_weighted_votes / pred_weight
      }else{
        rel_other = 0
      }
      #rel_other = ifelse(rel_other==0,0,log(rel_other)) 
      neighbor_info <- data.frame(
        other_score = rel_other,
        neighbor = paste(neighbors,collapse=","),
        distance = paste(round(distances, 3),collapse=","),
        #weight = paste(round(weights, 3),collapse=","),
        label = paste(neighbor_labels,collapse=","),
        weighted_votes = paste(weighted_votes,collapse=","),
        neighbors_other = neighbors_other,
        other_weighted_votes = other_weighted_votes,
        total_w = total_weight,
        pred_w = pred_weight
        
      )
      all_neighbors = bind_rows(all_neighbors, neighbor_info)
    }
    
  }
  to_return = data.frame(
    predicted_label = preds,
    confidence = confs)
  if(track_neighbors){
    #print(dim(to_return))
    #print(dim(all_neighbors))
    #check for any missing points
    if(!nrow(to_return)==nrow(all_neighbors)){
      print("mismatch in row number for to_return and all_neighbors")
      print(nrow(to_return))
      print(nrow(all_neighbors))
      print(head(to_return))
      print(head(all_neighbors))
      stop("")
    }
    to_return = bind_cols(to_return,all_neighbors)
  }
  return(to_return)
}


