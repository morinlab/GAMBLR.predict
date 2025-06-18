

#' Assemble genetic features for UMAP input
#'
#' This function assembles a matrix of genetic features for each sample, including mutation status,
#' aSHM counts, and structural variant status for BCL2, BCL6, and MYC. It supports both genome and capture sequencing types.
#'
#' @param these_samples_metadata Data frame with sample metadata, must include seq_type and sample_id.
#' @param metadata_columns Columns in metadata to use for SV status (default: c("bcl2_ba","bcl6_ba","myc_ba")).
#' @param genes Vector of gene symbols to include.
#' @param synon_genes Vector of gene symbols for synonymous mutations.
#' @param maf_with_synon MAF data frame including synonymous mutations.
#' @param hotspot_genes Vector of hotspot genes.
#' @param sv_value Value to assign for SV presence (default: 3).
#' @param synon_value Value to assign for synonymous mutations (default: 1).
#' @param coding_value Value to assign for coding mutations (default: 2).
#'
#' @return Matrix of assembled features for each sample.
#' @export
assemble_genetic_features <- function(these_samples_metadata,
                metadata_columns = c("bcl2_ba","bcl6_ba","myc_ba"),
                genes,
                synon_genes,
                maf_with_synon,
                hotspot_genes,
                genome_build = "grch37",
                sv_value = 3,
                synon_value = 1,
                coding_value = 2,
                include_ashm = TRUE){
  if(include_ashm){
      #TODO: ensure this supports both genome builds correctly
    some_regions = GAMBLR.utils::create_bed_data(
                  GAMBLR.data::grch37_ashm_regions,
                  fix_names = "concat",
                  concat_cols = c("gene","region"),
                  sep="-")
    #make the aSHM count matrix and combine if necessary
    if ("genome" %in% these_samples_metadata$seq_type){
    ashm_matrix_genome <- get_ashm_count_matrix(
    regions_bed = some_regions,
    this_seq_type = "genome",
    these_samples_metadata = these_samples_metadata
    )

    colnames(ashm_matrix_genome) = gsub("-.+","",colnames(ashm_matrix_genome))
    }
    if ("capture" %in% these_samples_metadata$seq_type){
    ashm_matrix_cap <- get_ashm_count_matrix(
    regions_bed = some_regions,
    this_seq_type = "capture",
    these_samples_metadata = these_samples_metadata
    )
    colnames(ashm_matrix_cap) = gsub("-.+","",colnames(ashm_matrix_cap))
    }
    if ("genome"  %in% these_samples_metadata$seq_type && "capture" %in% these_samples_metadata$seq_type){
      ashm_matrix = bind_rows(ashm_matrix_genome, ashm_matrix_cap)
    }else if("genome" %in% these_samples_metadata$seq_type){
      ashm_matrix = ashm_matrix_genome
    }else if("capture" %in% these_samples_metadata$seq_type){
      ashm_matrix = ashm_matrix_cap
    }else{
      stop("no eligible seq_type provided in these_samples_metadata")
    }
}
  
  
  status_with_silent = get_coding_ssm_status(
    these_samples_metadata = these_samples_metadata,
    maf_data = maf_with_synon,
    include_hotspots = TRUE,
    genes_of_interest = "MYD88",
    include_silent_genes = synon_genes[synon_genes %in% genes],
    gene_symbols = genes
  ) 
  status_with_silent = status_with_silent %>% column_to_rownames("sample_id")

  status_without_silent = get_coding_ssm_status(
    these_samples_metadata = these_samples_metadata,
    maf_data = maf_with_synon,
    include_hotspots = TRUE,
    genes_of_interest = "MYD88",
    include_silent = FALSE,
    gene_symbols = genes
  ) 

  status_without_silent = status_without_silent %>% column_to_rownames("sample_id")

  if(any(! colnames(status_with_silent) %in% colnames(status_without_silent))){
    print(colnames(status_with_silent)[!colnames(status_with_silent) %in% colnames(status_without_silent)])
    stop("some columns are missing from the status_without_silent matrix")
  }
  if(include_ashm){
    ashm_matrix = select(ashm_matrix, any_of(colnames(status_with_silent))) %>% select(any_of(synon_genes))

    #fill in gaps from aSHM (other non-coding variants in the genes)

    missing = status_with_silent[rownames(ashm_matrix),
                 colnames(ashm_matrix)]==0 & 
    ashm_matrix[rownames(ashm_matrix),
        colnames(ashm_matrix)] == synon_value

    status_with_silent[rownames(missing),
           colnames(missing)] = synon_value

  }
    

  status_combined = status_with_silent + status_without_silent
  if(coding_value == 1){
    status_combined[status_combined > 1] = 1
  }
  #TODO: generalize this to use the column names provided in metadata_columns
  bcl2_id = these_samples_metadata[which(these_samples_metadata$bcl2_ba=="POS"),] %>% pull(sample_id)
  bcl6_id = these_samples_metadata[which(these_samples_metadata$bcl6_ba=="POS"),] %>% pull(sample_id)
  myc_id = these_samples_metadata[which(these_samples_metadata$myc_ba=="POS"),] %>% pull(sample_id)

  status_combined[,"BCL2_SV"] = 0
  status_combined[bcl2_id,"BCL2_SV"] = sv_value

  status_combined[,"MYC_SV"] = 0
  status_combined[myc_id,"MYC_SV"] = sv_value

  status_combined[,"BCL6_SV"] = 0
  status_combined[bcl6_id,"BCL6_SV"] = sv_value
  return(status_combined)

}

old_assemble_genetic_features <- function(these_samples_metadata,
                              metadata_columns = c("bcl2_ba","bcl6_ba","myc_ba"),
                              genes,
                              synon_genes,
                              maf_with_synon,
                              hotspot_genes,
                              sv_value = 3,
                              synon_value = 1,
                              coding_value = 2){
  #TODO: ensure this supports both genome builds correctly
  some_regions = GAMBLR.utils::create_bed_data(
                              GAMBLR.data::grch37_ashm_regions,
                              fix_names = "concat",
                              concat_cols = c("gene","region"),
                              sep="-")
  #make the aSHM count matrix and combine if necessary
  if ("genome" %in% these_samples_metadata$seq_type){
    ashm_matrix_genome <- get_ashm_count_matrix(
     regions_bed = some_regions,
     this_seq_type = "genome",
     these_samples_metadata = these_samples_metadata
    )

    colnames(ashm_matrix_genome) = gsub("-.+","",colnames(ashm_matrix_genome))
  }
  if ("capture" %in% these_samples_metadata$seq_type){
    ashm_matrix_cap <- get_ashm_count_matrix(
     regions_bed = some_regions,
     this_seq_type = "capture",
     these_samples_metadata = these_samples_metadata
    )
    colnames(ashm_matrix_cap) = gsub("-.+","",colnames(ashm_matrix_cap))
  }
  if ("genome"  %in% these_samples_metadata$seq_type && "capture" %in% these_samples_metadata$seq_type){
    ashm_matrix = bind_rows(ashm_matrix_genome, ashm_matrix_cap)
  }else if("genome" %in% these_samples_metadata$seq_type){
    ashm_matrix = ashm_matrix_genome
  }else if("capture" %in% these_samples_metadata$seq_type){
    ashm_matrix = ashm_matrix_cap
  }else{
    stop("no eligible seq_type provided in these_samples_metadata")
  }
  
  status_with_silent = get_coding_ssm_status(
    these_samples_metadata = these_samples_metadata,
    maf_data = maf_with_synon,
    include_hotspots = TRUE,
    genes_of_interest = "MYD88",
    include_silent_genes = synon_genes,
    gene_symbols = genes
  ) 
  status_with_silent = status_with_silent %>% column_to_rownames("sample_id")

  status_without_silent = get_coding_ssm_status(
    these_samples_metadata = these_samples_metadata,
    maf_data = maf_with_synon,
    include_hotspots = TRUE,
    genes_of_interest = "MYD88",
    include_silent = FALSE,
    gene_symbols = genes
  ) 
  status_without_silent = status_without_silent %>% column_to_rownames("sample_id")

  status_combined = status_with_silent + status_without_silent
  if(coding_value == 1){
    status_combined[status_combined > 1] = 1
  }

  #fill in gaps from aSHM (other non-coding variants in the genes)

  missing = status_with_silent[rownames(ashm_matrix),
                               colnames(ashm_matrix)]==0 & 
    ashm_matrix[rownames(ashm_matrix),
                colnames(ashm_matrix)] == synon_value

  status_with_silent[rownames(missing),
                     colnames(missing)] = synon_value
  
  bcl2_id = these_samples_metadata[which(these_samples_metadata$bcl2_ba=="POS"),] %>% pull(sample_id)
  bcl6_id = these_samples_metadata[which(these_samples_metadata$bcl6_ba=="POS"),] %>% pull(sample_id)
  myc_id = these_samples_metadata[which(these_samples_metadata$myc_ba=="POS"),] %>% pull(sample_id)

  status_combined[,"BCL2_SV"] = 0
  status_combined[bcl2_id,"BCL2_SV"] = sv_value

  status_combined[,"MYC_SV"] = 0
  status_combined[myc_id,"MYC_SV"] = sv_value

  status_combined[,"BCL6_SV"] = 0
  status_combined[bcl6_id,"BCL6_SV"] = sv_value
  return(status_combined)

}

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
                              target_weight=0.5,
                              calc_dispersion = FALSE,
                              algorithm = "tumap"){
  
  # Function to compute mean (or median) pairwise distance within a group
  pairwise_dispersion <- function(df_group) {
    coords <- as.matrix(df_group[, c("V1", "V2")])
    dists <- dist(coords)  # Euclidean pairwise distances
    return(median(dists))    
  }
  if("sample_id" %in% colnames(df)){
    df = df %>% column_to_rownames(var = "sample_id")
  }
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
  if (!missing(metadata)) {
    df_sample_ids <- rownames(df)
    metadata_sample_ids <- metadata[[join_column]]
  
    # overlaps and mismatches
    shared_ids <- intersect(df_sample_ids, metadata_sample_ids)
    missing_from_metadata <- setdiff(df_sample_ids, metadata_sample_ids)
    missing_from_df <- setdiff(metadata_sample_ids, df_sample_ids)
  
    # shared rows
    df <- df[shared_ids, , drop = FALSE]
    metadata <- metadata %>% filter((!!sym(join_column)) %in% shared_ids)
  
    # otherwise message
    message(length(missing_from_metadata), " samples were missing from metadata.")
    message(length(missing_from_df), " samples in metadata were dropped as missing from the feature matrix.")
    message("using ", length(shared_ids), " samples provided in metadata and having features.")
  }
  
  if(missing(umap_out)){
    if(missing(target_column)){
      if(algorithm == "umap"){
        umap_out = umap2(df %>% as.matrix(),
                         n_neighbors = n_neighbors,
                         min_dist = min_dist,
                         metric = metric,
                         ret_model = ret_model,
                         n_epochs=n_epochs,
                         a=1.8956,
                         b = 0.806,
                         approx_pow=TRUE,
                         init=init,
                         seed = seed,
                         n_threads = 1, #IMPORTANT: n_threads must not be changed because it will break reproducibility
                         batch = TRUE,
                         n_sgd_threads = 1,
                         rng_type = "deterministic") # possibly add rng_type = "deterministic"  
        
      }else if(algorithm == "tumap"){
        umap_out = tumap(df %>% as.matrix(),
                         n_neighbors = n_neighbors,
                         metric = metric,
                         ret_model = ret_model,
                         n_epochs=n_epochs,
                         init=init,
                         seed = seed,
                         n_threads = 1,
                         batch = TRUE,
                         n_sgd_threads = 1,
                         rng_type = "deterministic")
      }else{
        stop("unsupported algorithm option")
      }
    }else{
      #supervised
      if(missing(metadata)){
        stop("metadata must be provided for supervised UMAP")
      }
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
                       n_threads = 1, #IMPORTANT: n_threads must not be changed because it will break reproducibility
                       y = metadata[[target_column]],
                       target_metric = target_metric,
                       target_weight = target_weight
                       ) # possibly add rng_type = "deterministic"
    }
    
  }else{
    umap_out = umap_transform(X=df,
                              model=umap_out$model,
                              seed=seed,
                              batch = TRUE,
                              n_threads = 1,
                              n_sgd_threads = 1)
    
  }
  
  if(!is.null(names(umap_out))){
    umap_df = as.data.frame(umap_out$embedding) %>% rownames_to_column(join_column)
  }else{
    umap_df = as.data.frame(umap_out) %>% rownames_to_column(join_column)
  }
  if(!missing(metadata)){
    umap_df = left_join(umap_df,metadata,by=join_column)
  }
  results = list()
  if(calc_dispersion){
      print("calculating pairwise dispersion")
      dispersion_df <- umap_df %>%
      group_by(lymphgen) %>%
      summarise(
        n = n(),
        mean_pairwise_distance = pairwise_dispersion(cur_data())
      ) %>%
      arrange(mean_pairwise_distance)
      results[["dispersion"]] = dispersion_df
  }
 
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
#' @param eval_group If desired, use this to specify which rows will be
#' evaluated and held out from training rather than using all samples. 
#' NOTE: this parameter will probably become deprecated!
#' @param umap_out The output of a previous run of make_and_annotate_umap.
#' If provided, the function will use this model to project the data
#' instead of re-running UMAP.
#' @param min_k Starting k for knn (Default: 3)
#' @param max_k Ending k for knn (Default: 33)
#' @param optimize_for_other Set to TRUE to optimize the threshold for 
#' classifying samples as "Other" based on the relative proportion of 
#' samples near the sample in UMAP space with the "Other" label. Rather than
#' treating Other as just another class, this will optimize the threshold for
#' a separate score that considers how many Other and non-Other samples are
#' in the neighborhood of the sample in question. This parameter will NOT change
#' the value in predicted_label. Instead, the predicted_label_optimized column
#' will contain the optimized label. Default: FALSE
#' 
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
#'
#' # Aim to maximize classification of samples into non-Other class
#' lymphgen_lyseq_no_other =  
#' GAMBLR.predict::DLBCLone_optimize_params(  
#'  dlbcl_status_combined_lyseq,
#'  dlbcl_meta_lyseq_train,
#'  min_k = 5,max_k=23,
#'  optimize_for_other = F,
#'  truth_classes = c("MCD","EZB","ST2","BN2"))
#'
#' # Aim to maximize balanced accuracy while allowing samples to be
#' # unclassified (assigned to "Other")
#' 
#' DLBCLone_lymphgen_lyseq_prime_opt =  
#' GAMBLR.predict::DLBCLone_optimize_params(  
#'  dlbcl_status_combined_lyseq_prime,
#'  dlbcl_meta_lyseq_train,
#'  min_k = 5,max_k=23,
#'  optimize_for_other = T,
#'  truth_classes = c("MCD","EZB","ST2","BN2","Other"))
#'

DLBCLone_optimize_params = function(combined_mutation_status_df,
                           metadata_df,
                           umap_out,
                           truth_classes = c("EZB",
                                             "MCD",
                                             "ST2",
                                             "N1",
                                             "BN2",
                                             "Other"),
                           optimize_for_other = FALSE,
                           eval_group = NULL,
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
      if(is.null(eval_group)){
        test_coords = filter(outs$df,lymphgen %in% truth_classes) %>% select(V1,V2)
        train_coords = filter(outs$df,lymphgen %in% truth_classes) %>% select(V1,V2)
        train_labels = filter(outs$df,lymphgen %in% truth_classes) %>% pull(lymphgen)
      }else{
        test_coords = filter(outs$df,dataset == eval_group) %>% select(V1,V2)
        train_coords = filter(outs$df,dataset != eval_group,lymphgen %in% truth_classes) %>% select(V1,V2)
        train_labels = filter(outs$df,dataset != eval_group,lymphgen %in% truth_classes) %>% pull(lymphgen)
      }

      pred_all = weighted_knn_predict_with_conf(
        train_coords = train_coords,
        train_labels = train_labels,
        test_coords = test_coords,
        k=ks,
        conf_threshold = -1, # Will threshold later
        na_label="Other",
        use_weights = use_w,
        ignore_top = ignore_top,
        verbose = verbose,
        track_neighbors = optimize_for_other
      )

      if(!"Other" %in% truth_classes){
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
          pred_other_all = weighted_knn_predict_with_conf(
            train_coords = train_coords,
            train_labels = train_labels,
            test_coords = test_coords,
            k=ks,
            conf_threshold = -1, # Will threshold later
            na_label="Other",
            use_weights = use_w,
            ignore_top = ignore_top,
            verbose = verbose
          )
        }
      } 

      for(k in ks){
        message(paste("K:",k))
        for (threshold in threshs){

          if(verbose){
            print(paste("k:",k,"threshold:",threshold,"use_weights:",use_w,"na_option:",na_option))
          }

          k_index <- which(ks == k)
          pred <- pred_all[[k_index]] %>%
            mutate(predicted_label = ifelse(confidence < threshold, "Other", predicted_label))
          
          if(is.null(eval_group)){
            xx_d = bind_cols(filter(outs$df,lymphgen %in% truth_classes) ,pred)
            train_d = filter(outs$df,lymphgen %in% truth_classes)
          }else{
            xx_d = bind_cols(filter(outs$df,dataset == eval_group) ,pred)
            train_d = filter(outs$df,dataset != eval_group,lymphgen %in% truth_classes)
          }
          
          # Always set true label factor levels to truth_classes
          xx_d$lymphgen = factor(xx_d$lymphgen, levels = truth_classes)

          # Set predicted labels to allow "Other"
          pred_levels = union(truth_classes, "Other")
          xx_d$predicted_label = factor(xx_d$predicted_label, levels = pred_levels)

          if(!"Other" %in% truth_classes && n_other > 0){
            pred_other <- pred_other_all[[k_index]] %>% 
              mutate(predicted_label = ifelse(confidence < threshold, "Other", predicted_label))

            if(!is.null(eval_group)){
              xx_o = bind_cols(
                filter(outs$df, dataset == eval_group, lymphgen == "Other"),
                pred_other
              )
            } else {
              xx_o = bind_cols(
                filter(outs$df, lymphgen == "Other"),
                pred_other
              )
            }

            xx_o$lymphgen = factor("Other", levels = pred_levels)
            xx_o$predicted_label = factor(xx_o$predicted_label, levels = pred_levels)
            xx_d = bind_rows(xx_d, xx_o)
          }

          true_factor = xx_d$lymphgen
          pred_factor = xx_d$predicted_label

          if (exclude_other_for_accuracy) {
            keep_idx = which(true_factor != "Other")
            true_factor = true_factor[keep_idx]
            pred_factor = pred_factor[keep_idx]
          }



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

            optimized_accuracy_and_thresh = optimize_outgroup(
              pred_factor,
              true_factor,
              xx_d$other_score,
              all_classes = truth_classes,
              maximize = maximize,
              exclude_other_for_accuracy = exclude_other_for_accuracy
            )
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
                            N1_bacc = unname(bal_acc["Class: N1"]),
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
            track_neighbors = optimize_for_other
  )
  pred = pred[[1]]

  if(optimize_for_other){
    pred = mutate(pred,predicted_label_optimized = ifelse(other_score > best_params$threshold_outgroup,
                                                          "Other",
                                                          predicted_label))
  }else{
    pred = mutate(pred,predicted_label_optimized = predicted_label)
  }
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
weighted_knn_predict_with_conf <- function(
  train_coords,
  train_labels,
  test_coords,
  k, # a vector
  epsilon = 0.1,
  conf_threshold = NULL,
  na_label = "Other",
  verbose = FALSE,
  use_weights = TRUE,
  ignore_top = FALSE,
  track_neighbors = TRUE,
  separate_other = TRUE #big change here. Other is considered separately for optimization
) { 
  if (nrow(train_coords)==0 || nrow(test_coords) == 0) {
    print("train_coords:")
    print(nrow(train_coords))
    print("test:")
    print(nrow(test_coords))
    stop("train_coords and test_coords must be data frames with at least one row")
  }

  train_labels = as.character(train_labels)
  max_k <- max(k)
  nn <- get.knnx(train_coords, test_coords, max_k)

  results_list <- list()

  for (curr_k in k) {
    if (verbose) cat("Processing k =", curr_k, "\n")

    preds <- character(nrow(test_coords))
    names(preds) <- rownames(test_coords) # *IMPORTANT* preserve sample_id

    confs <- numeric(nrow(test_coords))
    names(confs) <- rownames(test_coords) # *IMPORTANT* preserve sample_id

    all_neighbors <- data.frame()

    for (i in 1:nrow(test_coords)){
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

      distances = distances + epsilon
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
      # num_other_neighbors = sum(neighbor_labels == "Other")
      other_dists = distances[neighbor_labels == "Other"]
      if(separate_other){
        valid = valid & neighbor_labels != "Other"
      }
            
      neighbor_labels <- neighbor_labels[valid]
      weights <- weights[valid]
      distances <- distances[valid]
      neighbors <- neighbors[valid]

      #now take the first k neighbors
      if(length(neighbor_labels) > curr_k){
        neighbor_labels = neighbor_labels[1:curr_k]
        weights = weights[1:curr_k]
        distances = distances[1:curr_k]
        neighbors = neighbors[1:curr_k]
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
      mean_other_dist = mean(others_distances)
      if (length(neighbor_labels) == 0) {
        preds[i] <- "Other"
        confs[i] <- 1
        if(track_neighbors){
          rel_other = 10
          neighbor_info <- data.frame(
            other_score = rel_other,
            neighbor_id = paste(rownames(train_coords)[neighbors],collapse=","),
            neighbor = paste(neighbors,collapse=","),
            distance = paste(round(distances, 3),collapse=","),
            label = paste(neighbor_labels,collapse=","),
            weighted_votes = "",
            neighbors_other = neighbors_other,
            other_weighted_votes = 0,
            mean_other_dist = mean_other_dist,
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
          neighbor_id = paste(rownames(train_coords)[neighbors],collapse=","),
          neighbor = paste(neighbors,collapse=","),
          distance = paste(round(distances, 3),collapse=","),
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
      confidence = confs
    )

    if(track_neighbors){
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
    results_list[[paste0("k_", curr_k)]] <- to_return 
  }
  return(results_list)
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
                                  label_column = "predicted_label_optimized"){

  circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
    r = diameter / 2
    tt <- seq(0,2*pi,length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
  }
  if(missing(training_predictions)){
    training_predictions = single_sample_prediction_output$anno_df
  }
  single_sample_prediction_output$prediction = filter(single_sample_prediction_output$prediction,sample_id==this_sample_id)
  #extract the sample_id for all the nearest neighbors with non-Other labels
  my_neighbours = filter(single_sample_prediction_output$prediction,
                         sample_id == this_sample_id) %>% 
                  pull(neighbor_id) %>% strsplit(.,",") %>% unlist()
  #set up links connecting each neighbor to the sample's point
  links_df = filter(training_predictions,sample_id %in% my_neighbours) %>% mutate(group=lymphgen)

  sample_row <- single_sample_prediction_output$anno_df %>%
    filter(sample_id == this_sample_id) %>%
    slice_head(n = 1) # ensures 1 sample_id selected. Duplicates may occur when ignore_top = TRUE.
  my_x <- sample_row$V1
  my_y <- sample_row$V2

  if(prediction_in_title){
    title = paste(this_sample_id,
                  pull(single_sample_prediction_output$prediction,
                       !!sym(label_column)))
    if(single_sample_prediction_output$prediction$predicted_label_optimized == "Other" && single_sample_prediction_output$prediction$predicted_label !="Other"){
      title = paste(title,"(",single_sample_prediction_output$prediction$predicted_label,")")
    }

  }else{
    title = this_sample_id
  }
  links_df = mutate(links_df,my_x=my_x,my_y=my_y)
  links_df = links_df %>% select(V1, V2, my_x, my_y, group) %>% mutate(length = sqrt((V1 - my_x)^2 + (V2 - my_y)^2))  # Euclidean distance
  
  
  pp=ggplot(mutate(training_predictions,group=lymphgen),
         aes(x=V1,y=V2,colour=group)) + 
    geom_point(alpha=0.8,size=0.5) + 
    geom_segment(data=links_df,aes(x=V1,y=V2,xend=my_x,yend=my_y),alpha=0.5)+
    scale_colour_manual(values=get_gambl_colours()) + 
    ggtitle(title)+
    theme_minimal()
  if(add_circle){
    #add a circle around the sample
    d = quantile(links_df$length, 1)*2.1  # Euclidean distances and adding a 10% spacer
    circle = circleFun(c(my_x,my_y),diameter=d,npoints=100)
    pp = pp + geom_path(data=circle,aes(x=x,y=y),colour="black",alpha=1,inherit.aes=FALSE)
  }
  return(pp)
}

#' Predict class for a single sample without using umap_transform and plot result of classification
#'
#' @param seed Random seed for reproducibility
#' @param test_df Data frame containing the mutation status of the test sample
#' @param train_df Data frame containing the mutation status of the training samples
#' @param train_metadata Metadata for training samples with truth labels in lymphgen column
#' @param umap_out UMAP output from a previous run. The function will use this model to project the data, useful
#' for reproducibility and for using the same UMAP model on different datasets.
#' @param best_params Data frame from DLBCLone_optimize_params with the best parameters
#' @param other_df Data frame containing the predictions for samples in the "Other" class
#' @param predict_training Set to TRUE to predict the projected training samples once stored as stored_train_prediction. Use 
#' stored_train_prediction for subsequent use wehn using the same training set. Set to to FALSE to predict training samples every run 
#' @param stored_train_prediction Data frame containing the projected train predictions from train samples
#' @param ignore_top Set to TRUE to avoid considering a nearest neighbor with
#' distance = 0. This is usually only relevant when re-classifying labeled
#' samples to estimate overall accuracy
#' @param truth_classes Vector of classes to use for training and testing. Default: c("EZB","MCD","ST2","N1","BN2")
#' @param drop_unlabeled_from_training Set to TRUE to drop unlabeled samples from the training data
#' @param make_plot Set to TRUE to plot the UMAP projection and predictions
#' @param annotate_accuracy Set to true to add labels with accuracy values
#' @param label_offset Length of the label offset for the accuracy labels
#' @param title1 additional argument
#' @param title2 additional argument
#' @param title3 additional argument
#'
#' @returns a list of data frames with the predictions, the UMAP input, the UMAP projected output, the model, and a ggplot object
#' @export
#'
#' @examples
#' predict_single_sample_DLBCLone(
#'    seed = 1234,
#'    test_df = test_df,
#'    train_df = train_df,
#'    train_metadata = train_metadata,
#'    umap_out = umap_out,
#'    best_params = best_params
#'    predictions_df = predictions_df,
#'    annotate_accuracy = TRUE
#' )
#'
predict_single_sample_DLBCLone <- function(
  test_df,
  train_df,
  train_metadata,
  umap_out,
  best_params,
  other_df,
  predict_training = FALSE, 
  stored_train_prediction = NULL, 
  ignore_top = FALSE,
  truth_classes = c("EZB","MCD","ST2","N1","BN2"),
  drop_unlabeled_from_training=TRUE,
  make_plot = TRUE,
  annotate_accuracy = FALSE,
  label_offset = 2,
  title1="GAMBL",
  title2="predicted_class_for_HighConf",
  title3 ="predicted_class_for_Other",
  seed = 12345
){

  if(nrow(test_df)>1){
    message("Warning: you have supplied more than one sample to test with. Will proceed with all")
  }

  trained_features = colnames(umap_out$features)

  train_df = train_df %>%
    column_to_rownames("sample_id") %>%
    select(all_of(trained_features)) %>% 
    rownames_to_column("sample_id")
  train_id <- train_df$sample_id

  test_df = test_df %>%
    column_to_rownames("sample_id") %>%
    select(all_of(trained_features)) %>%
    rownames_to_column("sample_id")
  test_id <- test_df$sample_id

  if(ignore_top){
    if(test_id %in% train_df$sample_id){
      combined_df = train_df 
      combined_metadata <- train_metadata
    }else{
      combined_df = bind_rows(train_df,test_df)

      # Make dummy metadata for test samples
      test_metadata <- data.frame(
        sample_id = test_df$sample_id,
        cohort = "Test-Sample",
        lymphgen = NA # or any placeholder?
      )
      # Ensure train_metadata has matching columns
      train_metadata <- train_metadata %>% 
        select(sample_id, cohort, lymphgen)
      combined_metadata <- bind_rows(train_metadata, test_metadata)
    }
  } else {
    # Drop overlaps to prevent rowname collisions
    dupes <- intersect(train_df$sample_id, test_df$sample_id)
    if(length(dupes) > 0){
      warning(paste("Removing", length(dupes),
                    "overlapping samples from train_df to avoid duplicated rownames.\n",
                    "Consider setting ignore_top = TRUE to avoid inaccurate high confidence, and to keep all training samples."))
      train_df <- train_df %>% filter(!sample_id %in% dupes)
      train_metadata <- train_metadata %>% filter(!sample_id %in% dupes)
    }
    combined_df <- bind_rows(train_df, test_df)

    # Make dummy metadata for test samples
    test_metadata <- data.frame(
      sample_id = test_df$sample_id,
      cohort = "Test-Sample",
      lymphgen = NA # or any placeholder?
    )
    # Ensure train_metadata has matching columns
    train_metadata <- train_metadata %>% 
      select(sample_id, cohort, lymphgen)
    combined_metadata <- bind_rows(train_metadata, test_metadata)
  }

  projection <- make_and_annotate_umap(
    df = combined_df,
    metadata = combined_metadata,
    umap_out = umap_out,
    ret_model = FALSE,
    seed = seed,
    join_column = "sample_id",
    na_vals = best_params$na_option
  )

  train_coords = dplyr::filter(
    projection$df,
    sample_id %in% train_id
  ) %>% 
    select(sample_id,V1,V2) %>%
    column_to_rownames("sample_id")

  train_labels = dplyr::filter(                    
    projection$df,
    sample_id %in% train_id
  ) %>% 
    select(sample_id,V1,V2) %>%
    left_join( 
      combined_metadata %>% select(
        sample_id, 
        lymphgen
      ), 
      by = "sample_id"
    ) %>%
    pull(lymphgen) 

  test_coords = dplyr::filter(
    projection$df,
    sample_id %in% test_id
  ) %>% 
    select(sample_id,V1,V2) %>%
    column_to_rownames("sample_id")

  if(predict_training && !is.null(stored_train_prediction)){
    train_pred = stored_train_prediction
  }else{
    train_pred = weighted_knn_predict_with_conf(
      train_coords = train_coords,
      train_labels = train_labels,
      test_coords = train_coords, # <- predicitng training on self
      k = best_params$k,
      conf_threshold = best_params$threshold,
      na_label = "Other",
      use_weights = best_params$use_w,
      ignore_top = ignore_top
    )
    train_pred <- as.data.frame(train_pred)
    prefix <- paste0("k_", best_params$k, ".")
    colnames(train_pred) <- sub(prefix, "", colnames(train_pred))
    train_pred = rownames_to_column(train_pred, var = "sample_id")

    if(best_params$threshold_outgroup > 0){ # optimize_for_other = T in DLBCLone_optimize_params
      train_pred = mutate(train_pred,predicted_label_optimized = ifelse(
        other_score > best_params$threshold_outgroup,
        "Other",
        predicted_label
      ))
    }else{ # optimize_for_other = F
      train_pred = mutate(train_pred,predicted_label_optimized = predicted_label)
    }   
  }
    
  test_pred = weighted_knn_predict_with_conf(
    train_coords = train_coords,
    train_labels = train_labels,
    test_coords = test_coords,
    k = best_params$k,
    conf_threshold = best_params$threshold,
    na_label = "Other",
    use_weights = best_params$use_w,
    ignore_top = ignore_top,
  )
  test_pred <- as.data.frame(test_pred)
  prefix <- paste0("k_", best_params$k, ".")
  colnames(test_pred) <- sub(prefix, "", colnames(test_pred))
  test_pred = rownames_to_column(test_pred, var = "sample_id")

  if(best_params$threshold_outgroup > 0){ # optimize_for_other = T in DLBCLone_optimize_params
    test_pred = mutate(test_pred,predicted_label_optimized = ifelse(
      other_score > best_params$threshold_outgroup,
      "Other",
      predicted_label
    ))
  }else{ # optimize_for_other = F
    test_pred = mutate(test_pred,predicted_label_optimized = predicted_label)
  }

  anno_umap = select(projection$df, sample_id, V1, V2)

  anno_out = left_join(test_pred,anno_umap,by="sample_id")

  anno_out = anno_out %>%
    mutate(
      label = as.character(paste(
        sample_id,
        predicted_label_optimized,
        round(confidence,3)
      )),
      V1 = as.numeric(V1),
      V2 = as.numeric(V2)
    )

  predictions_train_df = left_join(train_pred, projection$df, by = "sample_id") 
  predictions_test_df = left_join(test_pred, projection$df, by = "sample_id")
  predictions_df = bind_rows(predictions_train_df, predictions_test_df)

  if(make_plot){
    title = paste0(
      "N_class:", best_params$num_classes," N_feats:",best_params$num_features,
      " k=",best_params$k," threshold=",best_params$threshold," bacc=",round(best_params$accuracy,3)
    )
    
    if("BN2" %in% truth_classes){
      print(best_params)
      acc_df = data.frame(
        lymphgen = c(
          "N1",
          "BN2",
          "EZB",
          "MCD",
          "ST2",
          "Other",
          "A53"
        ),
        accuracy = c(
          best_params$N1_bacc,
          best_params$BN2_bacc,
          best_params$EZB_bacc,
          best_params$MCD_bacc,
          best_params$ST2_bacc,
          best_params$Other_bacc,
          best_params$A53_bacc
        )
      )
    }else if("C1" %in% truth_classes){
      acc_df = data.frame(
        lymphgen = c(
          "C1",
          "C2",
          "C3",
          "C4",
          "C5"
        ),
        accuracy = c(
          best_params$C1_bacc,
          best_params$C2_bacc,
          best_params$C3_bacc,
          best_params$C4_bacc,
          best_params$C5_bacc
        )
      )
    }else{
      stop("no labels to add?")
    }

    # Add the predicted labels for Other (unclassified) cases, if provided
    if(!missing(other_df)){
      in_df = bind_rows(
        mutate(predictions_df,dataset=title2,lymphgen=predicted_label_optimized),
        mutate(other_df,dataset=title3,lymphgen=predicted_label_optimized)
      )
    }else{
      in_df = bind_rows(
        mutate(predictions_df,dataset=title2,lymphgen=predicted_label_optimized)
      )
    }

    pp = ggplot(in_df) +
      geom_point(aes(x=V1,y=V2,colour=lymphgen),alpha=0.8) +
      scale_colour_manual(values=get_gambl_colours()) +
      facet_wrap(~dataset,ncol=1) +
      theme_Morons() + ggtitle(title)

    if(annotate_accuracy){
      #add labels and set nudge direction based on what quadrant each group sits in
      centroids = filter(predictions_df,predicted_label_optimized %in% truth_classes) %>%
        group_by(predicted_label_optimized) %>%
        summarise(mean_V1=median(V1),mean_V2=median(V2)) %>%
        mutate(nudge_x=sign(mean_V1),nudge_y = sign(mean_V2)) %>%
        mutate(lymphgen=predicted_label_optimized)
        centroids = left_join(centroids,acc_df) %>%
        mutate(label=paste(lymphgen,":",round(accuracy,3)))

      pp = pp + 
        geom_label_repel(
          data=filter(centroids,nudge_y < 0, nudge_x < 0),
          aes(x=mean_V1,y=mean_V2,label=label),
          fill="white",
          size=5,
          nudge_y = -1 * label_offset, 
          nudge_x = -1 * label_offset
        ) +
        geom_label_repel(
          data=filter(centroids,nudge_y < 0, nudge_x > 0),
          aes(x=mean_V1,y=mean_V2,label=label),
          size=5,
          nudge_y = -1 * label_offset, 
          nudge_x = 1 * label_offset
        ) +
        geom_label_repel(
          data=filter(centroids,nudge_y > 0, nudge_x < 0),
          aes(x=mean_V1,y=mean_V2,label=label),
          size=5,
          nudge_y = 1 * label_offset, 
          nudge_x = -1 * label_offset
        ) +
        geom_label_repel(
          data=filter(centroids,nudge_y > 0, nudge_x > 0),
          aes(x=mean_V1,y=mean_V2,label=label),
          fill="white",
          size=5,
          nudge_y = 1 * label_offset, 
          nudge_x = 1 * label_offset
        ) +
        geom_label_repel(
          data = anno_out,
          aes(x=V1,y=V2,label=label),
          nudge_y = 1 * label_offset, 
          nudge_x = 1 * label_offset,
          colour="red"
        ) 
    }
  }

  return(list(
    prediction = test_pred, 
    train_prediction = train_pred,
    umap_input = umap_out$features, 
    model=umap_out$model,
    plot = pp,
    anno_df = predictions_df,
    projection = projection$df 
  ))
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
                                 add_labels = FALSE){
  
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
  
  p = ggplot(df,
             aes(x=V1,y=V2,colour=!!sym(colour_by),label=cohort)) + geom_point(alpha=0.8) + 
    scale_colour_manual(values=cols) + theme_Morons() + 
    guides(colour = guide_legend(nrow = 1))
  if(add_labels){
    p = p + geom_label_repel(data=labels,aes(x=median_x,y=median_y,label=!!sym(colour_by)))
  }
  ggMarginal(p,groupColour = TRUE,groupFill=TRUE)
  
}
