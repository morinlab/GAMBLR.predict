
#' Optimize the threshold for classifying samples as "Other"
#'
#' Performs a post-hoc evaluation of the classification of a sample as one of
#' the main classes vs the outgroup/unclassified label "Other" and returns the
#' optimal threshold for classifying a sample as "Other" based on the ground
#' truth provided in the true_labels vector. It evaluates the performance
#' of the classifier using a range of thresholds and returns the best threshold
#' based on the specified metric (balanced accuracy or accuracy). 
#' 
#' NOTE: This function is not generally meant to be called directly but rather is
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
                            exclude_other_for_accuracy = FALSE,
                            cap_classification_rate = 1,
                            verbose = FALSE,
                            other_class = "Other"){
  rel_thresholds = seq(1,10,0.1)
  sens_df = data.frame()
  acc_df = data.frame()
  predictions = data.frame(predicted_label=as.character(predicted_labels),
                           true_label=as.character(true_labels))

  for(threshold in rel_thresholds){
      predictions_new = mutate(predictions,
                               predicted_label = ifelse(other_score < threshold,
                                                        predicted_label,
                                                        other_class))
      all_acc = report_accuracy(predictions_new,
        truth="true_label",
        pred="predicted_label",
        drop_other=FALSE)
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
      bal_acc$harmonic_mean = all_acc$harmonic_mean
      bal_acc$classification_rate = all_acc$classification_rate
      acc_df = bind_rows(acc_df,bal_acc)
      sn <- conf_matrix$byClass[, "Sensitivity"]  
      sn$average_sensitivity = mean(sn)
      sn$threshold = threshold
      sn$harmonic_mean = all_acc$harmonic_mean
      sn$classification_rate = all_acc$classification_rate  
      sens_df = bind_rows(sens_df,sn)
  }
  
  
  if(maximize %in% c("balanced_accuracy","accuracy")){
    acc_df = filter(acc_df, classification_rate <= cap_classification_rate)
    best = slice_head(arrange(acc_df,desc(average_accuracy)),n=1)
  }else if(maximize == "harmonic_mean"){
    acc_df = filter(sens_df, classification_rate <= cap_classification_rate)
    best = slice_head(arrange(sens_df,desc(harmonic_mean)),n=1)
  }else{
    best = slice_head(arrange(sens_df,desc(average_sensitivity)),n=1)

  }
  
  return(best)
}



#' Process KNN Vote Strings and Scores for Classification
#'
#' This function processes the raw neighbor label strings and weighted vote scores from k-nearest neighbor (KNN) classification results.
#' It computes per-class neighbor counts, weighted scores, and identifies the top group by count and score for each sample.
#' The function also supports custom logic for handling the "Other" class, including vote multipliers and purity requirements.
#' NOTE: This is a helper function and is not intended to be called directly by the user
#'
#' @param df Data frame containing kNN results, including columns with neighbor labels and weighted votes.
#' @param raw_col Name of the column containing the comma-separated neighbor labels (default: "label").
#' @param group_labels Character vector of all possible class labels to consider (default: c("EZB", "MCD", "ST2", "BN2", "N1", "Other")).
#' @param vote_labels_col Name of the column containing the comma-separated neighbor labels for weighted votes (default: "vote_labels").
#' @param k Number of neighbors used in kNN (required).
#' @param other_vote_multiplier Multiplier for the "Other" class when determining if a sample should be reclassified as "Other" (default: 2).
#' @param score_purity_requirement Minimum ratio of top group score to "Other" score to assign a sample to the top group (default: 1).
#' @param weighted_votes_col Name of the column containing the comma-separated weighted votes (default: "weighted_votes").
#'
#' @return Data frame with additional columns for per-class neighbor counts, scores, top group assignments, and summary statistics for each sample.
#'
#' @details
#' - Computes the number of neighbors for each class and the sum of weighted votes per class.
#' - Identifies the top group by count and by weighted score, and applies custom logic for the "Other" class if present.
#' - Adds columns for counts, scores, top group, top group score, score ratios, and optimized group assignments.
#' - Designed for downstream use in DLBCLone and similar kNN-based classification workflows.
#'
#' @examples
#' # Example usage:
#' # result <- process_votes(knn_output_df, k = 7)
#'
#' @export
process_votes <- function(df,
                          raw_col = "label",
                          group_labels = c("EZB", "MCD", "ST2", "BN2", "N1", "Other"),
                          vote_labels_col = "vote_labels",
                          weighted_votes_col = "weighted_votes",
                          #vote_labels_col = "label",
                          k,
                          other_vote_multiplier = 2,
                          score_purity_requirement = 1,
                          
                          other_class = "Other",
                          optimize_for_other = TRUE,
                          debug = FALSE) {  
  if(missing(k)){
    stop("k value is required")
  }
  score_thresh = 2 * k

  count_labels_in_string <- function(string, labels) {
    tokens <- stringr::str_split(string, ",")[[1]]
    map_int(labels, ~ sum(tokens == .x))
  }

  extract_weighted_scores <- function(label_str, vote_str, labels) {
    lbls  <- stringr::str_split(label_str, ",")[[1]]
    votes <- as.numeric(stringr::str_split(vote_str, ",")[[1]])
    map_dbl(labels, ~ sum(votes[lbls == .x])) %>%
      set_names(paste0(labels, "_score"))
  }

  get_top_score_group <- function(label_str, vote_str, labels) {
    lbls  <- stringr::str_split(label_str, ",")[[1]]
    votes <- as.numeric(stringr::str_split(vote_str, ",")[[1]])
    scores_by_label <- set_names(map_dbl(labels, ~ sum(votes[lbls == .x])), labels)
    top    <- names(scores_by_label)[which.max(scores_by_label)]
    value  <- scores_by_label[[top]]
    list(top_score_group = top, top_group_score = value)
  }

  df_out <- df %>%
    mutate(.id = row_number()) %>%
    rowwise() %>%
    mutate(
      counts = list(
        set_names(
          count_labels_in_string(.data[[raw_col]], group_labels),
          paste0(group_labels, "_NN_count")
        )
      ),
      top_group = {
        cnts <- count_labels_in_string(.data[[raw_col]], group_labels)
        group_labels[which.max(cnts)]
      },
      scores = list(
        if (!is.null(vote_labels_col) && !is.null(weighted_votes_col)) {
          extract_weighted_scores(
            .data[[vote_labels_col]],
            .data[[weighted_votes_col]],
            group_labels
          )
        } else {
          set_names(rep(0, length(group_labels)), paste0(group_labels, "_score"))
        }
      ),
      score_summary = list(
        if (!is.null(vote_labels_col) && !is.null(weighted_votes_col)) {
          get_top_score_group(
            .data[[vote_labels_col]],
            .data[[weighted_votes_col]],
            group_labels
          )
        } else {
          list(top_score_group = NA_character_, top_group_score = NA_real_)
        }
      )
    ) %>%
    ungroup()
    

    df_out = df_out %>%
    unnest_wider(counts) %>%
    unnest_wider(scores) %>%
    unnest_wider(score_summary)
    long_version = df_out #save for debugging if necessary
    
    df_out = df_out %>%
    rowwise() %>%
    mutate(
      top_group_count = get(paste0(top_group, "_NN_count"))
    ) %>%
    ungroup()

  if (!is.null(other_class) && other_class %in% group_labels) {  
    if ("neighbors_other" %in% colnames(df)) {
      df_out <- df_out %>%
        mutate(!!sym(paste0(other_class, "_count")) := neighbors_other)
    }
    if ("other_weighted_votes" %in% colnames(df)) {
      df_out <- df_out %>%
        mutate(!!sym(paste0(other_class, "_score")) := other_weighted_votes)
    }
    if (all(c("top_group_count", paste0(other_class, "_count")) %in% colnames(df_out))) {
      df_out <- df_out %>%
        mutate(by_vote = top_group_count) %>%
        mutate(by_vote_opt = ifelse(top_group_count * other_vote_multiplier > !!sym(paste0(other_class, "_count")), top_group, other_class))
    }
    df_out <- mutate(df_out,
                     by_score = top_score_group,
                     score_ratio = top_group_score / !!sym(paste0(other_class, "_score")),
                     by_score_opt = ifelse(score_ratio > score_purity_requirement | top_group_score > score_thresh, top_score_group, other_class))
  } else {
    # Fallback: if no "other" class exists, keep by_score/by_vote as top_group
    df_out <- mutate(df_out,
                     by_score = top_score_group,
                     score_ratio = NA_real_,
                     by_score_opt = top_score_group,
                     by_vote_opt = top_group)
  }
  if(debug){
    return(list(processed=df_out, long_version = long_version))
  }else{
    return(df_out)
  }
  
}



#' Optimize Purity Threshold for Classification Assignment
#'
#' This function searches for the optimal purity threshold to assign samples to their predicted class or to "Other" based on the score ratio in processed kNN vote results.
#' It iteratively tests a range of purity thresholds, updating the predicted class if the score ratio meets or exceeds the threshold, and computes the accuracy for each threshold.
#' The function returns the best accuracy achieved and the corresponding purity threshold.
#' NOTE: This is a helper function and is not intended to be called directly by the user
#'
#' @param processed_votes Data frame output from `process_votes`, containing at least the columns for score ratio, by_score_opt, and the relevant prediction and truth columns.
#' @param prediction_column Name of the column in `processed_votes` to update with the optimized prediction.
#' @param truth_column Name of the column in `processed_votes` containing the true class labels.
#'
#' @return A list with two elements: `best_accuracy` (numeric, the highest accuracy achieved) and `best_purity_threshold` (numeric, the threshold at which this accuracy was achieved).
#'
#' @details
#' - For each threshold in the range 0.1 to 0.95 (step 0.05), the function updates the prediction column to assign the class from `by_score_opt` if the score ratio meets the threshold, otherwise assigns "Other".
#' - Accuracy is computed as the proportion of correct assignments (diagonal of the confusion matrix).
#' - The function is intended for use in optimizing classification purity in kNN-based workflows, especially when distinguishing between confident class assignments and ambiguous ("Other") cases.
#'
#' @import caret
#' @examples
#' # Example usage:
#' # result <- optimize_purity(processed_votes, prediction_column = "pred_label", truth_column = "true_label")
#'
#' @export
optimize_purity <- function(optimized_model_object,
                            vote_df, 
                            mode, 
                            optimize_by = "balanced_accuracy", #allowed: harmonic_mean, overall_accuracy
                            truth_column, 
                            all_classes = c("MCD","EZB","BN2","N1","ST2","Other"),
                            k,
                            cap_classification_rate = 1,
                            exclude_other_for_accuracy = FALSE,
                            other_class = "Other",
                            optimize_for_other = TRUE) {  

  out_column = "DLBCLone_wo"

  if(!missing(optimized_model_object)){
    if(!is.null(optimized_model_object$best_params)){
      if(!missing(k)){
        message("k is provided in the optimized_model_object, ignoring the k parameter")
      }
      k = optimized_model_object$best_params$k
    }else{
      stop("optimized_model_object must contain best_params with k value")
    }
    vote_df = optimized_model_object$predictions
  }

  some_classes <- if (!is.null(other_class) && other_class %in% all_classes) {
    all_classes[all_classes != other_class]
  } else {
    all_classes
  }
  score_thresh = 2 * k

  processed_votes <- process_votes(vote_df,
                                   group_labels = all_classes,
                                   k = k,
                                   score_purity_requirement = 0.5,
                                   other_class = other_class,
                                   optimize_for_other = optimize_for_other)    
  if(!truth_column %in% colnames(processed_votes)){
    stop("truth_column must be a column in processed_votes")
  }
  best_accuracy <- 0
  best_purity_threshold <- 0

  processed_votes <- mutate(processed_votes, !!sym(truth_column) := as.character(!!sym(truth_column))) 
  
  if (!is.null(other_class) && other_class %in% all_classes) {
    no_other_df <- processed_votes %>%
      filter(!!sym(truth_column) != other_class)
  } else {
    no_other_df <- processed_votes
  }

  for(purity_threshold in seq(3, 0, -0.05)){
    updated_votes <- processed_votes %>%
      mutate(!!sym(out_column) := if (!is.null(other_class) && other_class %in% all_classes) {
        ifelse(score_ratio >= purity_threshold | top_group_score > score_thresh, by_score, other_class)
      } else {
        by_score
      })

    updated_no_other_df <- no_other_df %>%
      mutate(!!sym(out_column) := if (!is.null(other_class) && other_class %in% all_classes) {
        ifelse(score_ratio >= purity_threshold | top_group_score > score_thresh, by_score, other_class)
      } else {
        by_score
      })

    if(!exclude_other_for_accuracy || is.null(other_class) || !(other_class %in% all_classes)){
      xx <- select(updated_votes, sample_id, !!sym(truth_column), !!sym(out_column)) %>%
        mutate(match = !!sym(truth_column) == !!sym(out_column)) %>%
        group_by(match) %>%
        summarise(concordant = sum(match == TRUE), discordant = sum(match == FALSE), .groups = "drop") %>%
        summarise(all_conc = sum(concordant), dis = sum(discordant), total = all_conc + dis, percent = 100 * sum(concordant) / total)
      
      updated_votes <- mutate(updated_votes, !!sym(out_column) := factor(!!sym(out_column), levels = all_classes))
      updated_votes <- mutate(updated_votes, !!sym(truth_column) := factor(!!sym(truth_column), levels = all_classes))
      #confusion_matrix <- table(updated_votes[[truth_column]], updated_votes[[out_column]])
      conf_matrix <- confusionMatrix(updated_votes[[out_column]], updated_votes[[truth_column]])
    } else {
      xx <- select(updated_no_other_df, sample_id, !!sym(truth_column), !!sym(out_column)) %>%
        filter(!!sym(truth_column) != other_class) %>%
        mutate(match = !!sym(truth_column) == !!sym(out_column)) %>%
        group_by(match) %>%
        summarise(concordant = sum(match == TRUE), discordant = sum(match == FALSE), .groups = "drop") %>%
        summarise(all_conc = sum(concordant), dis = sum(discordant), total = all_conc + dis, percent = 100 * sum(concordant) / total)

      updated_no_other_df <- mutate(updated_no_other_df,
        !!sym(out_column) := factor(!!sym(out_column), levels = all_classes),
        !!sym(truth_column) := factor(!!sym(truth_column), levels = all_classes))
      confusion_matrix <- table(updated_no_other_df[[truth_column]], updated_no_other_df[[out_column]])
      conf_matrix <- confusionMatrix(updated_no_other_df[[out_column]], updated_no_other_df[[truth_column]])
    }

 
    acc_details <- report_accuracy(updated_votes,
      truth = truth_column,
      pred = out_column)
    
    bal_acc <- mean(conf_matrix$byClass[, "Balanced Accuracy"], na.rm = TRUE)
    classification_rate = acc_details$classification_rate
    if(classification_rate <= cap_classification_rate){    
      if(optimize_by == "harmonic_mean"){
        accuracy = acc_details$harmonic_mean
      }else if( optimize_by == "overall_accuracy"){
        accuracy = conf_matrix$overall[["Accuracy"]]

      }else{
        accuracy = bal_acc
      }


      if(accuracy > best_accuracy){
        best_accuracy <- accuracy
        best_purity_threshold <- purity_threshold
      }
    }else{
      message(paste0("skipping purity threshold ",purity_threshold," because classification rate ",round(classification_rate,3),
                     " exceeds cap of ",cap_classification_rate))
    }
  }

  best_out <- processed_votes %>%
    mutate(!!sym(out_column) := if (!is.null(other_class) && other_class %in% all_classes) {
      ifelse(score_ratio >= best_purity_threshold | top_group_score > score_thresh, by_score, other_class)
    } else {
      by_score
    })

  if(missing(optimized_model_object)){
    optimized_model_object <- list()
    optimized_model_object$best_params <- list(
      k = k,
      purity_threshold = best_purity_threshold,
      accuracy = best_accuracy,
      num_classes = length(unique(best_out$predicted_label)),
      num_features = ncol(best_out) - 3, 
      seed = 12345
    )  
  }
  optimized_model_object$predictions <- best_out
  optimized_model_object$best_accuracy <- best_accuracy  
  optimized_model_object$best_purity_threshold <- best_purity_threshold
  optimized_model_object$maximize <- optimize_by
  optimized_model_object$score_thresh <- score_thresh
  optimized_model_object$exclude_other_for_accuracy <- exclude_other_for_accuracy

  return(optimized_model_object)
}


#' Summarize SSM (Somatic Single Nucleotide Mutation) Status Across Samples
#'
#' This function summarizes the mutation status for a set of genes
#' across multiple samples, separating mutations by class for genes specified
#' by \code{separate_by_class_genes} and counting the number of hits
#' in each sample per mutation category.
#'
#' @param maf_df A data frame containing mutation annotation format (MAF) data,
#' with at least the following columns:
#'   \code{Hugo_Symbol}, \code{Variant_Classification}, and \code{Tumor_Sample_Barcode}.
#' @param silent_maf_df (Optional) A separate data frame containing silent mutation data if
#' the user doesn't want to pull silent mutation status from \code{maf_df}.
#' This argument is useful when you want to combine mutations from
#' the output of get_coding_ssm and get_ssm_by_region or get_ssm_by_gene
#' @param these_samples_metadata A data frame containing metadata for the samples,
#' with at least a \code{sample_id} column.
#' Any sample that does not have a matching sample_id in these_samples_metadata will be dropped. 
#' @param genes_of_interest A character vector of gene symbols to include
#' in the summary. If missing, defaults to all Tier 1 B-cell lymphoma genes.
#' @param synon_genes (Optional) A character vector of gene symbols for which
#' synonymous mutations should be included.
#' @param separate_by_class_genes (Optional) A character vector of
#' gene symbols for which mutations should be separated by class
#' (e.g., "Nonsense_Mutation", "Missense_Mutation").
#' @param count_hits Logical; if \code{TRUE}, counts the number of mutations
#' per gene per sample. If \code{FALSE} (default), only presence/absence is recorded.
#'
#' @return A wide-format data frame (matrix) with samples as rows and
#' mutation types as columns. Each cell contains either the count of
#' mutations (if \code{count_hits = TRUE}) or a binary indicator (0/1)
#' for mutation presence.
#'
#' @details
#' - Mutations are grouped and optionally separated by mutation
#' class for genes specified in \code{separate_by_class_genes} 
#' - Synonymous mutations can be counted as another separate feature for genes specified by \code{synon_genes} genes.
#' - The function simplifies mutation annotations and pivots the data to a wide format suitable for downstream analysis.
#'
#' @examples
#' # A basic example, using only the output of get_all_coding_ssm
#' # Since the only non-coding class this function handles is Silent,
#' # we will be missing most non-coding types such as Intron, UTR, Flank
#' \dontrun{
#' sample_metadata = get_sample_metadata() %>% filter(seq_type!= "mrna")
#' 
#' maf_data = get_all_coding_ssm(sample_metadata)
#' 
#' mutation_matrix <- summarize_all_ssm_status(
#'   maf_df = maf_data,
#'   these_samples_metadata = sample_metadata,
#'   genes_of_interest = c("TP53", "SGK1", "BCL2"),
#'   synon_genes = c("BCL2"),
#'   separate_by_class_genes = c("TP53","SGK1"),
#'   count_hits = FALSE
#' )
#' }
#'
#' @import dplyr tidyr tibble readr
#' @export
summarize_all_ssm_status <- function(maf_df,
                                     these_samples_metadata,
                                     genes_of_interest,
                                     synon_genes,
                                     silent_maf_df,
                                     separate_by_class_genes = NULL, 
                                     count_hits = FALSE){
  if(missing(genes_of_interest)){
    message("defaulting to all Tier 1 B-cell lymphoma genes")
    genes_of_interest = filter(GAMBLR.data::lymphoma_genes,
      DLBCL_Tier==1 | FL_Tier == 1 | BL_Tier == 1 ) %>% 
      pull(Gene) %>% unique()
  }
  
  maf_df = filter(maf_df,
    Hugo_Symbol %in% genes_of_interest)
  if(missing(silent_maf_df)){
    silent_maf_df = maf_df
  }else{
    silent_maf_df = filter(silent_maf_df,Hugo_Symbol %in% genes_of_interest)
  }
  if(!missing(these_samples_metadata)){
    maf_df = filter(maf_df,
      Tumor_Sample_Barcode %in% these_samples_metadata$sample_id)
    silent_maf_df = filter(silent_maf_df,
                           Tumor_Sample_Barcode %in% these_samples_metadata$sample_id)
  }
  if(!missing(synon_genes)){
    if(any(!synon_genes %in% silent_maf_df$Hugo_Symbol)){
      missing = synon_genes[!synon_genes %in% silent_maf_df$Hugo_Symbol]
      not_missing = synon_genes[synon_genes %in% silent_maf_df$Hugo_Symbol]
      
      message("Warning: Some synonymous genes have no mutations in silent_maf_df")
      message(paste(missing,collapse=", "))
    }
    maf_nonsilent = filter(maf_df,  Variant_Classification %in% vc_nonSynonymous)
    maf_silent = filter(silent_maf_df, ! Variant_Classification %in% vc_nonSynonymous, Hugo_Symbol %in% synon_genes)
    maf_df = bind_rows(maf_silent,maf_nonsilent)
  }
  #Simplify annotations
  silent_types = c("Silent","Intron","5'UTR","3'UTR","5'Flank","3'Flank")
  nonsense_types = c("Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins","Splice_Site")
  missense_types = c("In_Frame_Del",
                              "In_Frame_Ins",
                              "Missense_Mutation",
                              "Splice_Region")
  gene_mutations = mutate(maf_df,
                    mutation_type=case_when(
                    Variant_Classification %in% nonsense_types ~ "Nonsense_Mutation",
                    Variant_Classification %in%  missense_types ~ "Missense_Mutation",
                    Variant_Classification %in% silent_types ~ "Silent",
                    TRUE ~ Variant_Classification)
                    ) %>%
                            mutate(mutation=paste(Hugo_Symbol,mutation_type,sep=":"))
  separated_coding_maf = filter(gene_mutations,
                                Hugo_Symbol %in% genes_of_interest,
                                Hugo_Symbol %in% separate_by_class_genes,
                                mutation_type != "Silent")
  
  unseparated_coding_maf = filter(gene_mutations,
                                  Hugo_Symbol %in% genes_of_interest, !Hugo_Symbol %in% separate_by_class_genes, mutation_type != "Silent") %>%
    mutate(mutation = paste(Hugo_Symbol,"Coding",sep=":"))
  separated_silent_maf = filter(gene_mutations,Hugo_Symbol %in% synon_genes, mutation_type == "Silent") %>%
    mutate(mutation = paste(Hugo_Symbol,"Silent",sep=":"))
  gene_mutations = bind_rows(separated_coding_maf, unseparated_coding_maf,separated_silent_maf)
  
  mutation_distinct = select(gene_mutations,mutation,Tumor_Sample_Barcode) %>% 
    mutate(mutated = 1) %>% 
    group_by(Tumor_Sample_Barcode,mutation,mutated) %>%
    count() %>% ungroup()
  if(count_hits){
    mutation_distinct = mutation_distinct %>% select(-mutated) %>% rename(mutated=n) %>% distinct()
  } else{
    mutation_distinct = mutation_distinct %>% select(-n) %>% distinct()
  }
  mutation_wide = pivot_wider(mutation_distinct, names_from = "mutation",values_from = "mutated",values_fill = 0) %>% 
  column_to_rownames("Tumor_Sample_Barcode")
 return(mutation_wide)
}




#' Check matrix against missing features.
#'
#' Operate on the matrix supplied by user and check it for any features missing
#'      compared to the provided set.
#'
#' @param incoming_matrix The incoming matrix to be checked for any missing
#'      features.
#' @param feature_set The feature set.
#' @return matrix where column names correspond to all features from the
#'      provided set
#'
check_for_missing_features <- function(
    incoming_matrix,
    feature_set
){
    # Check if any features are missing
    missing_features <- setdiff(
        feature_set,
        colnames(incoming_matrix)
    )

    if(length(missing_features)>0){
        message(
            "ATTENTION: Not all features are available in the data!"
        )
        message(
            paste0(
                "A total of ",
                length(
                   missing_features
                ),
                " features are missing:"
            )
        )
        message(
            paste(
                missing_features,
                collapse=", "
            )
        )
        message(
            "They will be set to 0, which may affect model performance."
        )
        incoming_matrix[,missing_features] <- 0

    }

    return(incoming_matrix)
}


#' Classify DLBCLs according to genetic subgroups of Chapuy et al.
#'
#' Use the feature weights from NMF model to assemble the binary matrix and
#'      classify DLBCL tumors based on C0-C5 system of Chapuy et al
#'
#' @param these_samples_metadata The metadata data frame that contains sample_id
#'      column with ids for the samples to be classified.
#' @param maf_data The MAF data frame to be used for matrix assembling. At least
#'      must contain the first 45 columns of standard MAF format.
#' @param seg_data The SEG data frame to be used for matrix assembling. Must be
#'      of standard SEG formatting, for example, as returned by get_cn_segments.
#' @param sv_data The SV data frame to be used for matrix assembling. Must be of
#'      standard BEDPE formatting, for example, as returned by get_combined_sv.
#' @param projection The projection of the samples. Only used to retrerive data
#'      through GAMBLR.data when it is not provided. Defaults to grch37.
#' @param output The output to be returned after prediction is done. Can be one
#'      of predictions, matrix, or both. Defaults to both.
#'
#' @return data frame, binary matrix, or both
#' @import dplyr readr GAMBLR.data
#'
classify_dlbcl_chapuy <- function(
    these_samples_metadata,
    maf_data,
    seg_data,
    sv_data,
    projection = "grch37",
    output = "both"
){
    # Assembling the feature matrix based on the guidance
    # non-synonymous mutations, 2; synonymous mutations, 1; no-mutation, 0;
    # high-grade CN gain [CN ≥ 3.7 copies], 2; low-grade CN gain [3.7 copies ≥ CN ≥ 2.2 copies], 1;
    # CN neutral, 0;
    # low-grade CN loss [1.1 ≤ CN ≤1.6 copies], 1; high-grade CN loss [CN ≤ 1.1 copies], 2;
    # chromosomal rearrangement present, 3; chromosomal rearrangement absent, 0
    chapuy_feature_matrix <- list()

    # Mutations matrix
    chapuy_feature_matrix$ssm_matrix <- maf_data %>%
        dplyr::filter(
            Hugo_Symbol %in% chapuy_features$ssm_features
        ) %>%
        dplyr::filter(
            Variant_Classification %in% c(
                "Silent",
                GAMBLR.data:::coding_class
            )
        ) %>%
        dplyr::select(
            Tumor_Sample_Barcode,
            Hugo_Symbol,
            Variant_Classification
        ) %>%
        dplyr::mutate(
            mutated = ifelse(
                Variant_Classification == "Silent",
                1,
                2
            )
        ) %>%
        dplyr::select(-Variant_Classification) %>%
        group_by(Tumor_Sample_Barcode,Hugo_Symbol) %>%
        dplyr::arrange(Tumor_Sample_Barcode, desc(mutated)) %>%
        dplyr::filter( # if both syn and nonsyn are present, prioritize nonsyn
            mutated==max(mutated)
        ) %>%
        group_by(Tumor_Sample_Barcode,Hugo_Symbol) %>%
        slice_head %>%
        ungroup %>%
        pivot_wider(
            names_from = "Hugo_Symbol",
            values_from = "mutated"
        ) %>%
        replace(is.na(.), 0) %>%
        column_to_rownames("Tumor_Sample_Barcode")

    chapuy_feature_matrix$ssm_matrix <- complete_missing_from_matrix(
        chapuy_feature_matrix$ssm_matrix,
        these_samples_metadata$sample_id
    )

    # CNV matrix
    if(projection=="grch37"){
        arm_coordinates <- GAMBLR.data::chromosome_arms_grch37
        cytoband_coordinates <- cytobands_grch37 %>%
            `names<-`(c("chr", "start", "end", "cytoband", "extra")) %>%
            dplyr::mutate(chr = gsub("chr", "", chr)) %>%
            dplyr::mutate(cytoband=paste0(chr,cytoband)) %>%
            dplyr::select(-extra)
    }else{
        arm_coordinates <- GAMBLR.data::chromosome_arms_hg38
        cytoband_coordinates <- cytobands_hg38 %>%
            `names<-`(c("chr", "start", "end", "cytoband", "extra")) %>%
            dplyr::mutate(cytoband=paste0(chr,cytoband)) %>%
            dplyr::select(-extra)
    }

    # First the arm features
    cnv_features_arm <- arm_coordinates %>%
        mutate(
            arm = paste0(
                chromosome,
                arm
            )
        ) %>%
        left_join(
            chapuy_features$cnv_features_arm,
            .,
            by="arm"
        )

    # Next, the cytoband features
    cnv_features_cytoband <- cytoband_coordinates %>%
        left_join(
            chapuy_features$cnv_features_cytoband,
            .,
            by="cytoband"
        )

    cnv_arms <- cool_overlaps(
            seg_data,
            cnv_features_arm,
            columns1 = c("chrom", "start", "end"),
            columns2 = c("chromosome", "start", "end")
        ) %>%
        dplyr::select(sample, arm, CNV, log.ratio) %>%
        dplyr::rename("feature"="arm")

    cnv_cytobands <-  cool_overlaps(
            seg_data,
            cnv_features_cytoband,
            columns1 = c("chrom", "start", "end"),
            columns2 = c("chr", "start", "end")
        ) %>%
        select(sample, cytoband, CNV, log.ratio) %>%
        rename("feature"="cytoband")

    chapuy_feature_matrix$cnv_matrix <- bind_rows(
            cnv_arms,
            cnv_cytobands
        ) %>%
        group_by(sample, feature, CNV) %>%
        summarise(
            featuremean = mean(log.ratio)
        ) %>%
        # get rid of neutrals
        dplyr::filter(
            !featuremean == 0
        ) %>%
        # ensure the same direction
        dplyr::filter(
            (featuremean>0 & CNV=="AMP") | (featuremean<0 & CNV=="DEL")
        ) %>%
        dplyr::mutate(
            CN = 2*2^featuremean,
            mutated = case_when(
                CNV=="AMP" & CN >=3.7 ~ 2,
                CNV=="AMP" & CN >=2.2 ~ 1,
                CN > 1.6 ~ 0,
                CNV=="DEL" & CN >1.1 ~ 1,
                CNV=="DEL" & CN <=1.1 ~ 2
                ),
            featurename = paste0(feature,":",CNV)
        ) %>%
        ungroup %>%
        dplyr::select(sample, mutated, featurename) %>%
        pivot_wider(
            names_from = "featurename",
            values_from = "mutated"
        ) %>%
        replace(is.na(.), 0) %>%
        column_to_rownames("sample")

    chapuy_feature_matrix$cnv_matrix <- complete_missing_from_matrix(
        chapuy_feature_matrix$cnv_matrix,
        these_samples_metadata$sample_id
    )

    # SV matrix
    chapuy_feature_matrix$sv_matrix <- sv_data %>%
        dplyr::filter(
            gene %in% chapuy_features$sv_features |
            partner %in% chapuy_features$sv_features
        ) %>%
        dplyr::mutate(
            feature = case_when(
                gene %in% chapuy_features$sv_features ~ paste0("SV:",gene),
                partner %in% chapuy_features$sv_features ~ paste0("SV:",partner)
            )
        ) %>%
        dplyr::mutate(
            mutated=3
        ) %>%
        distinct(
            tumour_sample_id, feature, mutated
        ) %>%
        ungroup %>%
        pivot_wider(
            names_from = "feature",
            values_from = "mutated"
        ) %>%
        replace(is.na(.), 0) %>%
        column_to_rownames("tumour_sample_id")

    chapuy_feature_matrix$sv_matrix <- complete_missing_from_matrix(
        chapuy_feature_matrix$sv_matrix,
        these_samples_metadata$sample_id
    )

    if("SV:CD274" %in% colnames(chapuy_feature_matrix$sv_matrix)){
        chapuy_feature_matrix$sv_matrix <- chapuy_feature_matrix$sv_matrix %>%
            dplyr::rename("SV:CD274/PDCD1LG2" = "SV:CD274")
    }


    # Generate complete matrix
    chapuy_feature_matrix$complete_matrix <- bind_cols(
        chapuy_feature_matrix$ssm_matrix,
        chapuy_feature_matrix$cnv_matrix,
        chapuy_feature_matrix$sv_matrix
    ) %>% as.data.frame

    # Check if any features are missing
    chapuy_feature_matrix$complete_matrix <- check_for_missing_features(
        chapuy_feature_matrix$complete_matrix,
        chapuy_features$feature_weights$Feature
    )

    # This is to ensure consistent ordering for a fool-proof downstream calculations
    chapuy_feature_matrix$complete_matrix <- chapuy_feature_matrix$complete_matrix %>%
        dplyr::select(chapuy_features$feature_weights$Feature)

    # If user only wants matrix, return it here and do not perform the
    # subsequent analysis
    if(output=="matrix"){
        return(chapuy_feature_matrix$complete_matrix)
    }

    # Classify the samples
    message("Assembled the matrix, classifying the samples ...")
    features_weights_matrix <- chapuy_features$feature_weights %>%
        column_to_rownames("Feature") %>%
        as.data.frame

    compute_cluster_probability <- function(Row) {
        ((Row %>% t) * features_weights_matrix) %>%
        colSums %>%
        as.data.frame %>%
        `names<-`(
            rownames(Row)
        )
    }

    predictions <- apply(
        chapuy_feature_matrix$complete_matrix,
        1,
        compute_cluster_probability
    )

    predictions <- do.call(
        cbind,
        predictions
    ) %>%
    as.data.frame %>%
    t %>% # The output is wide so convert it to have 1 row/sample
    as.data.frame

    predictions <- predictions %>%
        rownames_to_column("sample_id") %>%
        rowwise() %>%
        dplyr::mutate(
            C0 = ifelse(
                sum(C1:C5)==0,
                1,
                0
            ),
            .before = "C1"
        ) %>%
        ungroup %>%
        column_to_rownames("sample_id")

    # Scale sum of weights to arrive to vote confidence in [0,1] range
    predictions <- predictions %>%
        mutate(
            across(everything(), ~ . / pmax(rowSums(across(everything())), 1))
    )

    # Layer in which cluster the sample belongs to
    # by taking the highest sum of weights
    predictions <- predictions %>%
        mutate(
            predict = colnames(.)[max.col(across(everything()))]
        )

    predictions <- predictions %>%
        rownames_to_column("sample_id")

    # Account for C0 samples, which will have all weights calculated as 0
    predictions <- predictions %>%
        as.data.frame %>%
        dplyr::rename(
            "Chapuy_cluster"="predict"
        )

    if(output == "predictions"){
        return(predictions)
    }else if (output == "both") {
        return(
            list(
                matrix = chapuy_feature_matrix$complete_matrix,
                predictions = predictions
            )
        )
    }else{
        stop(
            paste0(
                "You requested to return ",
                output,
                ", which is not supported.\n",
                "Please specify one of matrix, predictions, or both."
            )
        )
    }

}


#' Classify DLBCLs according to genetic subgroups of Lacy et al.
#'
#' Use the random forest model to classify DLBCL tumors based on system of
#'      Lacy et al
#'
#' @param these_samples_metadata The metadata data frame that contains sample_id
#'      column with ids for the samples to be classified.
#' @param maf_data The MAF data frame to be used for matrix assembling. At least
#'      must contain the first 45 columns of standard MAF format.
#' @param seg_data The SEG data frame to be used for matrix assembling. Must be
#'      of standard SEG formatting, for example, as returned by get_cn_segments.
#' @param sv_data The SV data frame to be used for matrix assembling. Must be of
#'      standard BEDPE formatting, for example, as returned by get_combined_sv.
#' @param projection The projection of the samples. Used to annotate hotspot
#'      SSM mutations and retreive coordinates for shm features. Defaults to
#'      grch37.
#' @param output The output to be returned after prediction is done. Can be one
#'      of predictions, matrix, or both. Defaults to both.
#' @param include_N1 Whether to set samples with NOTCH1 truncating mutations to
#'      N1 group as described in Runge et al (2021). Defaults to FALSE.
#'
#' @return data frame, binary matrix, or both
#' @rawNamespace import(randomForest, except = c("combine"))
#' @import dplyr readr
#'
classify_dlbcl_lacy <- function(
    these_samples_metadata,
    maf_data,
    seg_data,
    sv_data,
    projection = "grch37",
    output = "both",
    include_N1 = FALSE
){

    # Assembling matrix of Lacy features
    lacy_feature_matrix <- list()

    maf_data <- annotate_hotspots(
        maf_data
    )

    # Mutations matrix
    lacy_feature_matrix$ssm <-
        tabulate_ssm_status(
            gene_symbols = lacy_features$ssm,
            these_samples_metadata = these_samples_metadata,
            maf_data = maf_data,
            include_hotspots = FALSE,
            genome_build = projection
        ) %>%
        column_to_rownames(
            "sample_id"
        )

    lacy_feature_matrix$ssm <- complete_missing_from_matrix(
        lacy_feature_matrix$ssm,
        these_samples_metadata$sample_id
    )

    # Hotspots matrix
    lacy_feature_matrix$hotspots <- maf_data %>%
        dplyr::filter(
            Hugo_Symbol %in% lacy_features$hotspots
        ) %>%
        dplyr::select(
            Tumor_Sample_Barcode,
            Hugo_Symbol,
            hot_spot
        ) %>%
        # when there is both hotspot and noncanonical, prioritize hotspots
        group_by(
            Tumor_Sample_Barcode,
            Hugo_Symbol
        ) %>%
        slice_head() %>%
        ungroup %>%
        # separate hotspot and noncanonical into different features
        mutate(
            Hugo_Symbol = ifelse(
                is.na(hot_spot),
                paste0(
                    Hugo_Symbol,
                    "_noncan"
                ),
                Hugo_Symbol
            ),
            mutate = 1
        ) %>%
        select(-hot_spot) %>%
        # convert to matrix format
        pivot_wider(
            names_from = "Hugo_Symbol",
            values_from = "mutate"
        ) %>%
        column_to_rownames(
            "Tumor_Sample_Barcode"
        ) %>%
        replace(is.na(.), 0) %>%
        # assure consistent ordering
        select(
            contains(
                "_noncan"
            ),
            everything()
        ) %>%
        rename(
            POU2F2_239 = any_of("POU2F2"),
            MYD88_265  = any_of("MYD88"),
            EZH2_646   = any_of("EZH2")
        )

    lacy_feature_matrix$hotspots <- complete_missing_from_matrix(
        lacy_feature_matrix$hotspots,
        these_samples_metadata$sample_id
    )

    # aSHM matrix
    ashm_features <- lacy_features[[paste0(projection, "_shm")]] %>%
        as.data.frame %>%
        mutate(
            start = feature_start,
            end = feature_end,
            name = Feature
        ) %>%
        select(
            chromosome,
            start,
            end,
            name
        )

    lacy_feature_matrix$shm <- cool_overlaps(
        maf_data,
        ashm_features,
        columns2 = colnames(ashm_features)[1:3]
    ) %>%
        group_by(Tumor_Sample_Barcode, name) %>%
        summarize(n = n()) %>%
        pivot_wider(
            id_cols = Tumor_Sample_Barcode,
            names_from = name,
            values_from = n,
            values_fill = 0
        ) %>%
        ungroup %>%
        distinct %>%
        column_to_rownames("Tumor_Sample_Barcode")

    lacy_feature_matrix$shm <- complete_missing_from_matrix(
        lacy_feature_matrix$shm,
        these_samples_metadata$sample_id
    )

    lacy_feature_matrix$shm[lacy_feature_matrix$shm>0] = 1

    # Amplifications were classed as driver events if they targeted a known
    # oncogene and resulted in predicted copy number of ≥ 6.
    # Deletions were classed as driver events if they targeted a known tumour
    # suppressor gene and resulted in heterozygous or homozygous loss
    lacy_feature_matrix$cnv <-
    base::get(paste0(projection, "_gene_coordinates")) %>%
        dplyr::filter(
            gene_name %in% lacy_features$cnv$Gene
        ) %>%
        left_join(
            lacy_features$cnv,
            .,
            by=c("Gene"="gene_name")
        ) %>%
        cool_overlaps(
            seg_data,
            .,
            columns1 = c("chrom", "start", "end"),
            columns2 = c("chromosome", "start", "end")
        ) %>%
        dplyr::mutate(
            CN = 2*2^log.ratio,
            mutated = 1
        ) %>%
        # drop neutrals
        dplyr::filter(
            !(CN > 1 & CN < 6)
        ) %>%
        # ensure the same direction
        dplyr::filter(
            (log.ratio>0 & CNV=="AMP") | (log.ratio<0 & CNV=="DEL")
        ) %>%
        dplyr::select(
            sample,
            Feature,
            mutated
        ) %>%
        # deduplicate
        group_by(
            sample,
            Feature
        ) %>%
        slice_head() %>%
        ungroup %>%
        # convert to matrix format
        pivot_wider(
            names_from = "Feature",
            values_from = "mutated"
        ) %>%
        column_to_rownames(
            "sample"
        ) %>%
        replace(is.na(.), 0)

    lacy_feature_matrix$cnv <- complete_missing_from_matrix(
        lacy_feature_matrix$cnv,
        these_samples_metadata$sample_id
    )

    # Generate complete matrix
    lacy_feature_matrix$complete <- bind_cols(
        lacy_feature_matrix$ssm %>% select(-any_of(lacy_features$hotspots)),
        lacy_feature_matrix$cnv,
        lacy_feature_matrix$shm,
        lacy_feature_matrix$hotspots
    ) %>% as.data.frame

    # Check if any features are missing
    lacy_feature_matrix$complete <- check_for_missing_features(
      lacy_feature_matrix$complete,
      lacy_features$all$Feature
    )

    # Handle the features where mutation or CNV is considered
    deduplicate_features <- function(Row, DataFrame) {
        DataFrame %>%
        dplyr::select(
            contains(
                Row[["Gene"]]
            )
        ) %>%
        dplyr::mutate(
            max_value = do.call(
                pmax, (.)
            )
        ) %>%
        dplyr::select(3) %>%
        `colnames<-`(Row[["Feature"]])
    }

    dedupliacted_ssm_cnv <- apply(
        lacy_features$cnv %>% dplyr::filter(Dual=="TRUE"),
        1,
        deduplicate_features,
        lacy_feature_matrix$complete
    )

    dedupliacted_ssm_cnv <- do.call(
        bind_cols,
        dedupliacted_ssm_cnv
    ) %>%
    as.data.frame

    # Now replace the original mut/CNV columns with augmented data
    lacy_feature_matrix$complete <-
        lacy_feature_matrix$complete %>%
        dplyr::select(-(
            lacy_features$cnv %>%
            dplyr::filter(Dual=="TRUE") %>%
            select(Gene,Feature) %>%
            as.list %>% unlist %>% unname
        )) %>%
        bind_cols(
            .,
            dedupliacted_ssm_cnv
        ) %>%
        as.data.frame %>%
        # This is to ensure consistent ordering for a fool-proof downstream calculations
        select(lacy_features$all$Feature)

    # If user only wants matrix, return it here and do not perform the
    # subsequent analysis
    if(output=="matrix"){
      return(lacy_feature_matrix$complete)
    }

    rf_matrix <- lacy_feature_matrix$complete

    # This is needed because RF does not handle colnames with special characters
    colnames(rf_matrix) <- (gsub(
        '_|-|\\.',
        '',
        colnames(rf_matrix)
        )
    )

    predictions = predict(
        RFmodel_Lacy,
        newdata=rf_matrix,
        type="Vote"
    ) %>%
    as.data.frame %>%
    rownames_to_column(
      "sample_id"
    )

    predictions = predict(
        RFmodel_Lacy,
        newdata=rf_matrix
    ) %>%
    as.data.frame %>%
    rownames_to_column(
      "sample_id"
    ) %>%
    dplyr::rename(
      "Lacy_cluster"="."
    ) %>%
    left_join(
      predictions,
      .,
      by="sample_id"
    ) %>%
    mutate(
      Lacy_cluster = as.character(Lacy_cluster)
    )

    # Manually set samples to N1 if user chooses to do so
    if(include_N1){
        n1_samples = maf_data %>%
            dplyr::filter(
                Hugo_Symbol=="NOTCH1"
            ) %>%
            dplyr::filter(
                grepl(
                "Frame_Shift",
                Variant_Classification
                )
            )  %>%
            pull(
                Tumor_Sample_Barcode
            )

        predictions <- predictions %>%
            dplyr::mutate(
                Lacy_cluster = ifelse(
                    sample_id %in% n1_samples,
                    "NOTCH1",
                    as.character(Lacy_cluster)
                )
            ) %>%
            dplyr::rename(
                "hmrn_cluster" = "Lacy_cluster"
            ) %>%
            dplyr::mutate(
                NOTCH1 = ifelse(
                    sample_id %in% n1_samples,
                    1,
                    0
                ),
                .before = "hmrn_cluster"
            )  %>%
            mutate(
                across(
                    where(is.numeric) & !c(NOTCH1),
                    ~ ifelse(NOTCH1 == 1, 0, .)
                )
            )
        if (output == "both") {
            return(
                list(
                    matrix = lacy_feature_matrix$complete,
                    predictions = predictions
                )
            )
        }

    }

    if(output == "predictions"){
      return(predictions)
    }else if (output == "both") {
      return(
        list(
          matrix = lacy_feature_matrix$complete,
          predictions = predictions
        )
      )
    }else{
      stop(
        paste0(
          "You requested to return ",
          output,
          ", which is not supported.\n",
          "Please specify one of matrix, predictions, or both."
        )
      )
    }

}


#' Construct LymphGenerator matrix
#'
#' Use the LymphGenerator features to construct matrix. Optionally, flatten some features having the same
#' biological information.
#'
#' @param these_samples_metadata The metadata data frame that contains sample_id column with ids for the samples to be classified. Must also contain column pathology.
#' @param maf_data The MAF data frame to be used for matrix assembling. At least must contain the first 45 columns of standard MAF format.
#' @param sv_data The SV data frame to be used for matrix assembling. Must be of standard BEDPE formatting, for example, as returned by get_combined_sv.
#' @param seg_data The SEG data frame to be used for matrix assembling. Must be of standard SEG formatting, for example, as returned by get_sample_cn_segments. Must be already adjusted for ploidy.
#' @param seq_type String of the seq type for the sample set.
#' @param projection String of projection of the samples. Only used to retrieve data through GAMBLR.data when it is not provided. Defaults to grch37.
#' @param output The output to be returned. Currently only matrix is supported.
#' @param drop_after_flattening Boolean on whether to remove features (rows) after flattening. Defaults to FALSE.
#' @return binary matrix
#' @import GAMBLR.data dplyr readr tibble GAMBLR.helpers
#'
classify_dlbcl_lymphgenerator <- function(
	these_samples_metadata,
    maf_data,
    sv_data,
    seg_data,
	seq_type = "genome",
	projection = "grch37",
    output = "matrix",
    drop_after_flattening = FALSE
){
    # Before computation check that user requested only matrix
    if(!output == "matrix"){
      stop(
        "At this point LymphGenerator mode can only support matrix generation"
      )
    }

    # Standardize the genome build
    projection <- handle_genome_build(projection)

    # Initiate placeholders
    features <- list()
    matrix <- list()

    if(!missing(seg_data)){
        # Collect CNV features
        oncogenes_bed <- lymphgenerator_features$CNV %>%
            dplyr::filter(
                genome_build == projection
            ) %>%
            dplyr::select(chromosome, start, end, symbol, direction)

        # drop low-level events
        seg_data <- seg_data %>%
            dplyr::filter(
                sample %in% these_samples_metadata$Tumor_Sample_Barcode
            ) %>%
            dplyr::mutate(
                CN = (2 * 2 ^ log.ratio)
            ) %>%
            dplyr::filter(
                abs(log.ratio) > 0.56 &
                !CN %in% c(1, 2, 3)
            )

        # Generate the CNV matrix
        matrix$cnv <- cool_overlaps(
            oncogenes_bed,
            seg_data,
            nomatch = TRUE,
            columns1 = c("chromosome", "start", "end"),
            columns2 = c("chrom", "start", "end")
        ) %>%
        dplyr::rename("Tumor_Sample_Barcode" = "sample") %>%
        dplyr::select(Tumor_Sample_Barcode, symbol, CN, direction) %>%
        dplyr::filter(
            (CN < 2 & direction == "LOSS") |
            (CN > 2 & direction %in% c("AMP", "GAIN"))
        ) %>%
        dplyr::mutate(
            mutated = 1,
            symbol = paste0(symbol, "_", direction)
        ) %>%
        dplyr::select(-c(CN, direction)) %>%
        dplyr::distinct(
            symbol,
            Tumor_Sample_Barcode,
            mutated
        ) %>%
        pivot_wider(
            .,
            names_from = "symbol",
            values_from = "mutated",
            values_fill = 0
        ) %>%
        column_to_rownames("Tumor_Sample_Barcode")

        matrix$cnv <- complete_missing_from_matrix(
            matrix$cnv,
            list_of_samples = these_samples_metadata$Tumor_Sample_Barcode
        )
    }

    if(!missing(sv_data)){
        # Collect SV features
        matrix$sv <- sv_data %>%
            dplyr::filter(
                tumour_sample_id %in% these_samples_metadata$Tumor_Sample_Barcode,
                gene %in% lymphgenerator_features$SV,
                !is.na(partner)
            ) %>%
            dplyr::mutate(
                feature = case_when(
                    gene == "BCL6" & partner == "MYC" ~ "MYC_BCL6",
                    gene == "MYC" & ! partner  == "IGH" ~ "MYC_SV",
                    gene == "MYC" ~ paste0(gene, "_", partner),
                    TRUE ~ paste0(gene, "_SV")
                    ),
                mutated = 1
            ) %>%
            dplyr::distinct(
                tumour_sample_id,
                feature,
                mutated
            ) %>%
            pivot_wider(
                .,
                names_from = "feature",
                values_from = "mutated",
                values_fill = 0
            ) %>%
            column_to_rownames("tumour_sample_id")

        matrix$sv <- complete_missing_from_matrix(
            matrix$sv,
            list_of_samples = these_samples_metadata$Tumor_Sample_Barcode
        )
    }

    # Collect HOTSPOT features
    matrix$hotspot <- annotate_hotspots(
        maf_data %>% filter(
            Hugo_Symbol %in% lymphgenerator_features$hotspot,
            Tumor_Sample_Barcode %in% these_samples_metadata$Tumor_Sample_Barcode),
        recurrence_min = 10
    )

    matrix$hotspot <- review_hotspots(matrix$hotspot) %>%
        filter(hot_spot == TRUE) %>%
        select(Tumor_Sample_Barcode, Hugo_Symbol, hot_spot)

    matrix$hotspot <- matrix$hotspot %>%
        mutate(hot_spot = 1) %>%
        distinct(
            Tumor_Sample_Barcode,
            Hugo_Symbol,
            .keep_all = TRUE
        ) %>%
        pivot_wider(
            .,
            names_from = "Hugo_Symbol",
            values_from = "hot_spot",
            values_fill = 0
        ) %>%
        column_to_rownames("Tumor_Sample_Barcode")

    matrix$hotspot <- complete_missing_from_matrix(
        matrix$hotspot,
        list_of_samples = these_samples_metadata$Tumor_Sample_Barcode
    )

    colnames(matrix$hotspot) <- paste0(
        colnames(matrix$hotspot),
        "HOTSPOT"
    )

    # Collect SSM features, including NFKBIZ 3'
    matrix$ssm <- tabulate_ssm_status(
        gene_symbols = lymphgenerator_features$SSM,
        these_samples_metadata = these_samples_metadata,
        maf_data = maf_data,
        include_hotspots = FALSE,
        include_silent = FALSE,
        genome_build = projection
    ) %>%
    column_to_rownames("sample_id")

    matrix$ssm <- complete_missing_from_matrix(
        matrix$ssm,
        list_of_samples = these_samples_metadata$Tumor_Sample_Barcode
    )

    if("NFKBIZ" %in% maf_data$Hugo_Symbol){
        matrix$nfkbiz <- maf_data %>%
            dplyr::filter(
                Hugo_Symbol == "NFKBIZ",
                Variant_Classification == "3'UTR",
                Tumor_Sample_Barcode %in% these_samples_metadata$Tumor_Sample_Barcode
        ) %>%
        dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol) %>%
        distinct() %>%
        dplyr::mutate(
            Hugo_Symbol = paste0(
                Hugo_Symbol,
                "_3UTR"
            ),
            mutated = 1
        ) %>%
        pivot_wider(
            .,
            names_from = "Hugo_Symbol",
            values_from = "mutated",
            values_fill = 0
        ) %>%
        dplyr::mutate(these_names = Tumor_Sample_Barcode) %>%
        column_to_rownames("these_names")

        matrix$nfkbiz <- complete_missing_from_matrix(
            matrix$nfkbiz,
            list_of_samples = these_samples_metadata$Tumor_Sample_Barcode
        ) %>%
        select(-Tumor_Sample_Barcode)

        matrix$ssm <- cbind(
            matrix$ssm,
            matrix$nfkbiz
        )
    }

    # Collect aSHM features
    message(
        paste(
            "Assembling the aSHM matrix for LymphGenerator using seq type",
            seq_type
        )
    )

    ashm_bed <- base::get(
        paste0(
            projection, "_ashm_regions"
        )
    ) %>%
    mutate(
        name = paste(
            gene, region, sep = "-"
        )
    ) %>%
    dplyr::filter(name %in% lymphgenerator_features$aSHM) %>%
    rowwise %>%
    mutate(
        chr_name = ifelse(
            projection == "grch37",
            gsub("chr", "", chr_name),
            chr_name
        )
    ) %>%
    ungroup %>%
    as.data.frame
    names(ashm_bed)[1:3] <- c("chrom", "start", "end")

    matrix$ashm <- cool_overlaps(
        maf_data,
        ashm_bed,
        columns2 = colnames(ashm_bed)[1:3]
    ) %>%
        group_by(Tumor_Sample_Barcode, name) %>%
        summarize(n = n()) %>%
        pivot_wider(
            id_cols = Tumor_Sample_Barcode,
            names_from = name,
            values_from = n,
            values_fill = 0
        ) %>%
        ungroup %>%
        distinct %>%
        column_to_rownames("Tumor_Sample_Barcode")


    matrix$ashm <- complete_missing_from_matrix(
        matrix$ashm,
        these_samples_metadata$sample_id
    )

    # Binarizing the matrix
    matrix$ashm <- matrix$ashm %>%
        rownames_to_column("Tumor_Sample_Barcode") %>%
        left_join(
            .,
            these_samples_metadata %>%
                select(Tumor_Sample_Barcode, pathology),
            by = "Tumor_Sample_Barcode"
        ) %>%
        column_to_rownames("Tumor_Sample_Barcode")

    # calculate average N of mutations in each aSHM site per pathology
    matrix$ashm_aggregated <- aggregate(
        . ~ pathology,
        data = matrix$ashm,
        FUN = mean
    ) %>%
    dplyr::relocate(
        pathology,
        .after = last_col()
    )

    # convert ashm counts and averages to long format
    matrix$ashm <- matrix$ashm %>%
        rownames_to_column("Tumor_Sample_Barcode") %>%
        pivot_longer(
            !c(Tumor_Sample_Barcode, pathology),
            names_to = "Region",
            values_to = "N_mut"
        )

    matrix$ashm_aggregated <- matrix$ashm_aggregated %>%
        pivot_longer(
            !pathology,
            names_to = "Region",
            values_to = "Average_mut"
        )

    # merge counts/averages together
    matrix$ashm <- left_join(
        matrix$ashm,
        matrix$ashm_aggregated,
        by = c("pathology", "Region")
    )

    matrix$ashm <- matrix$ashm %>%
        dplyr::mutate(
            Average_mut = ifelse(
                pathology == "DLBCL",
                Average_mut + 3,
                Average_mut + 1)
        ) %>%
        dplyr::mutate(
            Average_mut = ifelse(
                seq_type == "capture",
                0,
                Average_mut)
        ) %>%
        dplyr::mutate(
            Feature = ifelse(
                N_mut > Average_mut,
                1,
                0)
        ) %>%
        dplyr::select(
            Tumor_Sample_Barcode,
            Region,
            Feature
        ) %>%
        pivot_wider(
            names_from = "Region",
            values_from = "Feature"
        ) %>%
        dplyr::arrange(
            match(
                Tumor_Sample_Barcode,
                these_samples_metadata$Tumor_Sample_Barcode)
        ) %>%
        column_to_rownames("Tumor_Sample_Barcode")

    # Drop helper matrices here before merge
    matrix <- matrix[!names(matrix) %in% c("nfkbiz", "ashm_aggregated")]

    matrix$full <- do.call(
        bind_cols,
        matrix
    ) %>%
    as.data.frame %>%
    t

    #####
    # Now flatten the features that we want to be squished
    for(i in 1:length(GAMBLR.predict:::names)){
        matrix$full <- flatten_feature(
            GAMBLR.predict:::names[i],
            GAMBLR.predict:::features[[i]],
            matrix$full
        )
    }

    matrix$full <- massage_matrix_for_clustering(
        matrix$full,
        blacklisted_cnv_regex = "3UTR|SV|HOTSPOT|MYC|BCL2|TP53BP1|intronic",
        min_feature_percent = 0
    )

    if(drop_after_flattening){
        matrix$full <- matrix$full[
            !(rownames(matrix$full) %in% unlist(features))
            ,
        ]
    }

    matrix$full <- matrix$full %>%
        t() %>%
        as.data.frame()

    return(matrix$full)
}


#' Harmonize different flavors of genome builds.
#'
#' Will process different genome build flavors and return it in consistent formatting used throughout this package.
#'
#' @param incoming_genome_build The string specifying the genome build that is about to be harmonized.
#' @return string
#'
#'
handle_genome_build <- function(
  incoming_genome_build
){

    if (incoming_genome_build %in% hg19_build_flavours){
        this_genome_build = "grch37"
    }else if(incoming_genome_build %in% hg38_build_flavours){
        this_genome_build = "hg38"
    }else{
        stop(
        paste(
            "The specified genome build",
            incoming_genome_build,
            "is not currently supported."
            )
        )
    }

    return(this_genome_build)
}

#' Flatten feature
#'
#' Return matrix (features in rows and samples in columns) where several features are
#' squished (flattened) into a single one. The number and order of features to be flattened
#' does not matter. New row with the new name specified by user will be created, and will be positive
#' (value of 1) if any of the features to flatten are positive.
#'
#' @param new_name String of the new feature name.
#' @param features_to_flatten Vector of features (row names) to be flattened.
#' @param incoming_data Matrix or data frame of features.
#' @return matrix
#' @import dplyr tidyselect
#'

flatten_feature <- function(
    new_name,
    features_to_flatten,
    incoming_data
){
    incoming_data <- incoming_data %>%
        as.data.frame %>%
        t() %>%
        as.data.frame

    incoming_data <- check_for_missing_features(
        incoming_data,
        features_to_flatten
    )

    flattened_data <- incoming_data %>%
        dplyr::mutate(
            across(
                {{ features_to_flatten }}, ~replace(., . == 0, NA)
            )
        ) %>%
        dplyr::mutate(
            !!new_name := coalesce(!!! syms ( features_to_flatten ))
        ) %>%
        dplyr::select(
            sort(tidyselect::peek_vars())
        ) %>%
        t %>%
        replace(is.na(.), 0)

    return(flattened_data)
}


#' @title Get Coding SSM Status.
#'
#' @description Tabulate mutation status (SSM) for a set of genes.
#'
#' @details This function takes a data frame (in MAF-like format) and converts
#' it to a binary one-hot encoded matrix of mutation status for either a set of
#' user-specified genes (via gene_symbols) or, if no genes are provided, default
#' to all lymphoma genes. The default behaviour is to assign each gene/sample_id
#' combination as mutated only if there is a protein coding mutation for that
#' sample in the MAF but this can be configured to use synonymous variants in
#' some (via include_silent_genes) or all (via include_silent) genes.
#' This function also has other filtering and convenience parameters giving
#' the user full control of the return. For more information, refer to the
#' parameter descriptions and examples.
#' Currently only the grch37 genome build is supported for hotspot annotation
#' and review for this version of the function.
#'
#' @param gene_symbols A vector of gene symbols for which the mutation status
#'      will be tabulated. If not provided, lymphoma genes will be returned
#'      by default.
#' @param these_samples_metadata The metadata for samples of interest to be
#'      included in the returned matrix. Only the column "sample_id" is
#'      required. If not provided, the example metadata is used as default.
#' @param maf_data data frame in maf format. Must be in the grch37 projection.
#' @param include_hotspots Logical parameter indicating whether hotspots object
#'      should also be tabulated. Default is TRUE.
#' @param keep_multihit_hotspot Logical parameter indicating whether to keep the
#'      gene annotation as mutated when the gene has both hot spot and
#'      non-hotspot mutation. Default is FALSE. If set to TRUE, will report the
#'      number of non-hotspot mutations instead of tabulating for just mutation
#'      presence.
#' @param review_hotspots Logical parameter indicating whether hotspots object
#'      should be reviewed to include functionally relevant mutations or rare
#'      lymphoma-related genes. Default is TRUE.
#' @param genes_of_interest A vector of genes for hotspot review. Currently only
#'      FOXO1, MYD88, and CREBBP are supported.
#' @param genome_build Reference genome build for the coordinates in the MAF
#'      file. The default is inferred from maf_data.
#' @param include_silent Logical parameter indicating whether to include silent
#'      mutations into coding mutations. Default is FALSE.
#' @param include_silent_genes Optionally, provide a list of genes for which the
#'      Silent variants to be considered. If provided, the Silent variants for
#'      these genes will be included regardless of the include_silent argument.
#' @param ... Any other parameter. These parameters will be ignored.
#'
#' @return A data frame with tabulated mutation status.
#'
#' @import dplyr tidyr GAMBLR.helpers
#' @export
#'
#' @examples
#' coding_tabulated_df = tabulate_ssm_status(
#'  maf_data = get_coding_ssm(),
#'  gene_symbols = c("EZH2","KMT2D","CREBBP","MYC")
#' )
#'
#'
#'
#' #all lymphoma genes from bundled NHL gene list
#' coding_tabulated_df = tabulate_ssm_status()
#'
#' #this example will fail because hg38 is not supported by this function (yet)
#' coding_tabulated_df = tabulate_ssm_status(maf_data=
#'                         get_coding_ssm(projection = "hg38"))
#' # Error in tabulate_ssm_status(maf_data = get_coding_ssm(projection = "hg38")) :
#' # Currently only grch37 projection (hg19 genome build) is supported.
#'
tabulate_ssm_status = function(
        gene_symbols,
        these_samples_metadata,
        maf_data,
        include_hotspots = TRUE,
        keep_multihit_hotspot = FALSE,
        review_hotspots = TRUE,
        genes_of_interest = c("FOXO1", "MYD88", "CREBBP"),
        genome_build,
        include_silent = FALSE,
        include_silent_genes,
        ...
    ){
    if(missing(maf_data)){
        stop("maf_data is required")
    }
    if("maf_data" %in% class(maf_data)){
        maf_data <- as.data.frame(maf_data)
    }
    # check the projection
    if(!genome_build == "grch37"){
        stop(
            "Currently only grch37 projection (hg19 genome build) is supported."
        )
    }

    if(missing(gene_symbols)){
        message(
            "No gene_symbols provided, defaulting to all lymphoma genes."
        )
        gene_symbols <- GAMBLR.data::lymphoma_genes$Gene
    }

    if(!missing(include_silent_genes)){
        message(
            strwrap(
                prefix = " ",
                initial = "",
                "Output will include all genes specified in gene_symbols
                and include_silent_genes parameters."
            )
        )
        gene_symbols <- c(
            gene_symbols,
            include_silent_genes
        ) %>%
        unique()
    }

    if(missing(these_samples_metadata)){
        stop(
            "Please provide sample metadata."
        )
    }

    coding_var <- c(
        "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del",
        "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation",
        "Nonstop_Mutation", "Splice_Region", "Splice_Site",
        "Targeted_Region", "Translation_Start_Site"
    )

    if(include_silent){
        message("Including Synonymous variants for all genes...")
        coding_var <- c(coding_var, "Silent")
    }

    if(missing(include_silent_genes)){
        coding_ssm <- maf_data %>%
            dplyr::filter(
                Variant_Classification %in% coding_var
            )
    } else {
        message(
            strwrap(
                prefix = " ",
                initial = "",
                "You have provided gene list with argument include_silent_genes.
                The Silent variants will be included even if the include_silent
                argument is set to FALSE.
                "
            )
        )
        coding_ssm <- maf_data %>%
            dplyr::filter(
                Variant_Classification %in% coding_var |
                (
                    Hugo_Symbol %in% include_silent_genes &
                    Variant_Classification == "Silent"
                )
            )
    }

    coding <- coding_ssm %>%
        dplyr::filter(
            Hugo_Symbol %in% gene_symbols
        ) %>%
        dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol) %>%
        dplyr::rename(
            "sample_id" = "Tumor_Sample_Barcode",
            "gene" = "Hugo_Symbol"
        ) %>%
        unique() %>%
        dplyr::mutate(mutated = 1)

    samples_table <- dplyr::select(
        these_samples_metadata,
        sample_id
    )
    wide_coding <- pivot_wider(
        coding,
        names_from = "gene",
        values_from = "mutated",
        values_fill = 0
    )
    all_tabulated <- left_join(
        samples_table,
        wide_coding
    )
    all_tabulated <- all_tabulated %>%
        replace(is.na(.), 0)

    # include hotspots if user chooses to do so
    if(include_hotspots){
        # first annotate
        annotated <- left_join(
            coding_ssm,
            GAMBLR.data::hotspots_annotations,
            by = c("Chromosome", "Start_Position")
        )

        # review for the supported genes
        if(review_hotspots){
            annotated = review_hotspots(
                annotated,
                genes_of_interest = genes_of_interest,
                genome_build = genome_build
            )
        }

        message("annotating hotspots")

        hotspots <- annotated %>%
            dplyr::filter(Hugo_Symbol %in% genes_of_interest) %>%
            dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol, hot_spot) %>%
            dplyr::rename(
                "sample_id" = "Tumor_Sample_Barcode",
                "gene" = "Hugo_Symbol"
            ) %>%
            dplyr::mutate(gene = paste0(gene, "HOTSPOT")) %>%
            unique() %>%
            dplyr::mutate(mutated = ifelse(hot_spot == "TRUE", 1, 0)) %>%
            dplyr::filter(mutated == 1) %>%
            dplyr::select(-hot_spot)

        # long to wide hotspots, samples are tabulated with 0 if no hotspot is detected
        wide_hotspots <- pivot_wider(
            hotspots,
            names_from = "gene",
            values_from = "mutated",
            values_fill = 0
        )
        # join with the ssm object
        all_tabulated <- left_join(
            all_tabulated,
            wide_hotspots
        )
        all_tabulated <- all_tabulated %>%
            replace(is.na(.), 0)

        all_tabulated <- all_tabulated %>%
            dplyr::select(where(~ any(. != 0)))

        all_tabulated <- as.data.frame(all_tabulated)
        # make SSM and hotspots non-redundant by giving priority
        # to hotspot feature and setting SSM to 0
        for (hotspot_site in colnames(wide_hotspots)[grepl("HOTSPOT", colnames(wide_hotspots))]){
            message(hotspot_site)
            this_gene <- gsub("HOTSPOT", "", hotspot_site)
            redundant_features <- all_tabulated %>%
                dplyr::select(starts_with(this_gene))

            # if not both the gene and the hotspot are present, go to
            # the next iteration
            if(ncol(redundant_features)!= 2) next
            message("OK")
            # if both gene and it's hotspot are in the matrix, give priority to hotspot feature
            all_tabulated[(all_tabulated[, this_gene] >0 & all_tabulated[, paste0(this_gene, "HOTSPOT")] == 1),][,c(this_gene, paste0(this_gene, "HOTSPOT"))][, this_gene] = 0

            # in case gene has both hotspot and another mutation in the same gene,
            # keep both tabulated as multihits
            if(keep_multihit_hotspot){
                # determine which samples have hot spot and another mutation in same gene
                multihits <- annotated %>%
                    dplyr::filter(Hugo_Symbol == this_gene) %>%
                    group_by(Tumor_Sample_Barcode) %>%
                    dplyr::mutate(n_mut = n()) %>%
                    dplyr::filter(
                        n_mut > 1
                    ) %>%
                    dplyr::distinct(Tumor_Sample_Barcode, n_mut, hot_spot) %>%
                    # account for cases with both hotspot and not hotspot to avoid
                    # double-counting the number of mutations
                    mutate_at(vars(hot_spot), ~replace_na(., "FALSE")) %>%
                    dplyr::mutate(
                        n_mut = ifelse(
                            hot_spot == "TRUE",
                            n_mut - 1,
                            n_mut
                        )
                    ) %>%
                    group_by(Tumor_Sample_Barcode) %>%
                    dplyr::arrange(n_mut) %>%
                    slice_head() %>%
                    ungroup %>%
                    select(-hot_spot)

                # Return the annotation of this gene to mutated in these samples
                all_tabulated <- all_tabulated %>%
                    left_join(
                        .,
                        multihits,
                        by = c("sample_id" = "Tumor_Sample_Barcode")
                    ) %>%
                    dplyr::mutate(
                        {{this_gene}} := ifelse(
                                !is.na(n_mut),
                                n_mut,
                                !!!syms(this_gene)
                            )
                    ) %>%
                    select(- n_mut)
            }

        }

    }
    all_tabulated <- distinct(all_tabulated)
    return(all_tabulated)

}
