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
#' @param verbose Defaults to FALSE
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
                include_ashm = TRUE,
                annotated_sv,
                include_GAMBL_sv= TRUE,
                review_hotspots = TRUE,
                verbose = FALSE){

  if(!"maf_data" %in% class(maf_with_synon)){
    warning(paste("maf_with_synon should be a maf_data object, but is not.",
                "Proceeding anyway with genome_build = ", genome_build))
  }else{
    genome_build = attr(maf_with_synon,"genome_build")
  }
  if(include_ashm){
    stopifnot( genome_build == "grch37", "other genome builds are not yet supported for aSHM")
    #TODO: consider shifting this function to GAMBLR.results
    #TODO: ensure this supports both genome builds correctly. This is currently hard-coded
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
  
  include_hotspots = ifelse(!missing(hotspot_genes), TRUE, FALSE)

  status_with_silent = get_coding_ssm_status(
    these_samples_metadata = these_samples_metadata,
    # drop all coding variants from this one
    maf_data = maf_with_synon,
    include_hotspots = include_hotspots,
    genes_of_interest = hotspot_genes,
    include_silent_genes = synon_genes[synon_genes %in% genes],
    gene_symbols = genes
  ) 
  status_with_silent = status_with_silent %>% column_to_rownames("sample_id")

  status_without_silent = get_coding_ssm_status(
    these_samples_metadata = these_samples_metadata,
    maf_data = maf_with_synon,
    include_hotspots = include_hotspots,
    genes_of_interest = hotspot_genes,
    include_silent = FALSE,
    gene_symbols = genes
  ) 

  status_without_silent = status_without_silent %>% 
    column_to_rownames("sample_id")

  if (any(! colnames(status_with_silent) %in% colnames(status_without_silent))){
    print(colnames(status_with_silent)[!colnames(status_with_silent) %in% colnames(status_without_silent)])
    stop("some columns are missing from the status_without_silent matrix")
  }
  # Instead of just relying on the MAF(s) supplied by the user
  if (include_ashm){
    ashm_matrix = select(ashm_matrix, any_of(colnames(status_with_silent))) %>% select(any_of(synon_genes))
    ashm_matrix[ashm_matrix>1]= 1
    
    if(verbose){
      print(head(colSums(ashm_matrix)))
      print(head(ashm_matrix[,c(1:10)]))
    }
    #fill in gaps from aSHM (other non-coding variants in the genes)

    missing = status_with_silent[rownames(ashm_matrix),
                 colnames(ashm_matrix)]==0 & 
    ashm_matrix[rownames(ashm_matrix),
        colnames(ashm_matrix)] > 0 
    
    
    fill = missing
    fill[]=0
    fill[missing] = synon_value
    status_with_silent[rownames(fill),
           colnames(fill)] = fill

  }

  #ensure all columns in status_with_silent are present in status_without_silent
  if (any(! colnames(status_without_silent) %in% colnames(status_with_silent))){
    print(colnames(status_without_silent)[!colnames(status_without_silent) %in% colnames(status_with_silent)])
    stop("some columns are missing from the status_with_silent matrix")
  }

  # ensure column order is identical in the two matrices

  status_with_silent = status_with_silent[,colnames(status_without_silent)]

  status_combined = status_with_silent + status_without_silent
  if(coding_value == 1){
    status_combined[status_combined > 1] = 1
  }
  #TODO: generalize this to use the column names provided in metadata_columns
  bcl2_id = these_samples_metadata[which(these_samples_metadata$bcl2_ba=="POS"),] %>% pull(sample_id)
  bcl6_id = these_samples_metadata[which(these_samples_metadata$bcl6_ba=="POS"),] %>% pull(sample_id)
  myc_id = these_samples_metadata[which(these_samples_metadata$myc_ba=="POS"),] %>% pull(sample_id)
 
  # include SV from provided annotated SV data frame
  if(!missing(annotated_sv) & include_GAMBL_sv){
    annotated_sv = filter(annotated_sv,tumour_sample_id %in% these_samples_metadata$sample_id)
    bcl2_sv_id = filter(annotated_sv,!is.na(partner),gene=="BCL2") %>% pull(tumour_sample_id)
    bcl6_sv_id = filter(annotated_sv,!is.na(partner),gene=="BCL6") %>% pull(tumour_sample_id)
    myc_sv_id = filter(annotated_sv,!is.na(partner),gene=="MYC") %>% pull(tumour_sample_id)
    n_b = length(bcl2_id)

    bcl2_id = unique(c(bcl2_id,bcl2_sv_id))
    n_b = length(bcl2_id)

    bcl6_id = unique(c(bcl6_id,bcl6_sv_id))
    myc_id = unique(c(myc_id,myc_sv_id))
  }
  

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
                            exclude_other_for_accuracy = FALSE,
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
#' @import uwot dplyr tidyr tibble readr
#' @export
#'
#' @examples
#' 
#' #library(GAMBLR.predict)
#'
#' all_full_status = readr::read_tsv(system.file("extdata/all_full_status.tsv",package = "GAMBLR.predict")) %>%
#'  tibble::column_to_rownames("sample_id")
#' dlbcl_meta = readr::read_tsv(system.file("extdata/dlbcl_meta_with_dlbclass.tsv",package = "GAMBLR.predict"))
#' my_umap <- make_and_annotate_umap(
#'   df=all_full_status,
#'   metadata=dlbcl_meta
#' )
#'
make_and_annotate_umap = function(df,
                              metadata,
                              umap_out,
                              truth_column = "lymphgen",
                              core_features = NULL,
                              core_feature_multiplier = 1.5,
                              hidden_features = NULL,
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
                              algorithm = "tumap",
                              make_plot = FALSE
                              ){
  
  if(!missing(umap_out)){
    model_provided = TRUE
  }else{
    model_provided = FALSE
  }
  if("sample_id" %in% colnames(df)){
    df = df %>% column_to_rownames(var = "sample_id")
  }
  original_n = nrow(df)
  if(na_vals == "to_zero"){
    df[is.na(df)] = 0
  }else if(na_vals == "drop"){
    if(any(sapply(df, is.factor))){
      numeric_cols = names(df)[!sapply(df, is.factor)]
      df <- df[, colSums(is.na(df[,numeric_cols])) == 0]
        rs = rowSums(df[,numeric_cols],na.rm=TRUE)
        dropped_rows = df[rs==0,]
        df = df[rs>0,]
    } else{
      df <- df[, colSums(is.na(df)) == 0]
        rs = rowSums(df,na.rm=TRUE)
        dropped_rows = df[rs==0,]
        df = df[rs>0,]
    }
    
  }
 

  if(missing(df)){
    stop("provide a data frame or matrix with one row for each sample and a numeric column for each mutation feature")
  }
  if(!is.null(core_features)){
    #print(paste(core_features,collapse=","))
    if(!is.numeric(core_feature_multiplier)){
      stop("core_feature_multiplier must be a numeric value")
    }
    if(!all(core_features %in% colnames(df))){
      stop("core_features must be a vector of column names in df")
    }
    #multiply the core features by the multiplier
    ncore = length(core_features)
    message(paste0("multiplying ",ncore," core features by ",core_feature_multiplier))
    df[core_features] = df[core_features] * core_feature_multiplier
  }
  if(!is.null(hidden_features)){
    if(!all(hidden_features %in% colnames(df))){
      stop("hidden_features must be a vector of column names in df")
    }

    message(paste0("dropping ",length(hidden_features)," hidden features"))
    df = df %>% select(-any_of(hidden_features))
  }
  no_feat_samples = NULL
  if(!missing(metadata)){
    if(nrow(df)< original_n){
      nrem = original_n-nrow(df)
      pct_rem = round(nrem/original_n*100,2)
      message(paste0("removed ",nrem," (",pct_rem,"%) rows from the data that had no features"))
      no_feat_samples = rownames(dropped_rows)
  
    }

    keep_rows = rownames(df)[rownames(df) %in% metadata[[join_column]]]
    df= df[keep_rows,]
    no_feat_metadata = filter(metadata,!!sym(join_column) %in% no_feat_samples)
    metadata= filter(metadata,!!sym(join_column) %in% rownames(df))
    
    message(paste("kept",nrow(metadata),"rows of the data that have features and match the metadata provided"))
  }
  
  if(missing(umap_out)){
    if(missing(target_column)){
      if(algorithm == "umap"){
        message("running umap2")
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
                         n_threads = 1,
                         batch = TRUE,
                         n_sgd_threads = 1,
                         rng_type = "deterministic") # possibly add rng_type = "deterministic"
        #IMPORTANT: n_threads must never be changed because it will break reproducibility  
        
      }else if(algorithm == "tumap"){
        
        message("running tumap")
        X = df
        umap_args = list(X=X,
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
      
        umap_out = tumap(X,
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
        message("done!")
      }else{
        stop("unsupported algorithm option")
      }
    }else{
      #supervised
      if(missing(metadata)){
        stop("metadata must be provided for supervised UMAP")
      }
      metadata[[target_column]] = factor(metadata[[target_column]])
      
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
                       target_weight = target_weight,
                       rng_type = "deterministic"
                       ) 
      #IMPORTANT: n_threads must never be changed because it will break reproducibility
    }
    

  }else{
    message("transforming each data point individually using the provided UMAP model. This will take some time.")
    umap_df = data.frame()
    for(sample in rownames(df)){
      this_row = df[sample,]
      this_umap_df = umap_transform(X=this_row,
                              model=umap_out$model,
                              seed=seed,
                              batch = TRUE,
                              n_threads = 1,
                              n_sgd_threads = 1)
      this_umap_df = as.data.frame(this_umap_df) 
      
      umap_df = bind_rows(umap_df,this_umap_df)
    }
    message("done")

    
  }
  
  
  if(model_provided){
    # model was generated here
    #message("model given to function")
    
    umap_df = as.data.frame(umap_df) %>% rownames_to_column(var=join_column)
    
  }else{
    umap_df = as.data.frame(umap_out$embedding) %>% rownames_to_column(join_column)
  }
  if(!missing(metadata)){
    umap_df = left_join(umap_df,metadata,by=join_column)
  }
  results = list()


  
  results[["df"]]=umap_df
  results[["features"]] = df
  results[["dropped_rows"]] = no_feat_samples
  if(!missing(metadata)){
    results[["total_samples_available"]] = nrow(metadata) + nrow(no_feat_metadata)
    results[["sample_metadata_no_features"]] = no_feat_metadata
  }
  
  if(ret_model){
    results[["model"]]= umap_out
  }
  if(make_plot){
    ms = make_umap_scatterplot(  umap_df,
                                  colour_by = truth_column,
                                  title = "UMAP projection")
    print(ms)
    results[["plot"]] = ms
  }
  return(results)
}

#' Process KNN Vote Strings and Scores for Classification
#'
#' This function processes the raw neighbor label strings and weighted vote scores from k-nearest neighbor (KNN) classification results.
#' It computes per-class neighbor counts, weighted scores, and identifies the top group by count and score for each sample.
#' The function also supports custom logic for handling the "Other" class, including vote multipliers and purity requirements.
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
                          k,
                          other_vote_multiplier = 2,
                          score_purity_requirement = 1,
                          weighted_votes_col = "weighted_votes",
                          other_class = "Other") {   # <--- NEW ARG
  if(missing(k)){
    stop("k value is required")
  }
  score_thresh = 2 * k

  count_labels_in_string <- function(string, labels) {
    tokens <- str_split(string, ",")[[1]]
    map_int(labels, ~ sum(tokens == .x))
  }

  extract_weighted_scores <- function(label_str, vote_str, labels) {
    lbls  <- str_split(label_str, ",")[[1]]
    votes <- as.numeric(str_split(vote_str, ",")[[1]])
    map_dbl(labels, ~ sum(votes[lbls == .x])) %>%
      set_names(paste0(labels, "_score"))
  }

  get_top_score_group <- function(label_str, vote_str, labels) {
    lbls  <- str_split(label_str, ",")[[1]]
    votes <- as.numeric(str_split(vote_str, ",")[[1]])
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
    ungroup() %>%
    unnest_wider(counts) %>%
    unnest_wider(scores) %>%
    unnest_wider(score_summary) %>%
    rowwise() %>%
    mutate(
      top_group_count = get(paste0(top_group, "_NN_count"))
    ) %>%
    ungroup()

  # Optional adjustments based on external columns (if present)
  if (!is.null(other_class) && other_class %in% group_labels) {   # <--- CONDITIONAL
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

  return(df_out)
}



#' Optimize Purity Threshold for Classification Assignment
#'
#' This function searches for the optimal purity threshold to assign samples to their predicted class or to "Other" based on the score ratio in processed kNN vote results.
#' It iteratively tests a range of purity thresholds, updating the predicted class if the score ratio meets or exceeds the threshold, and computes the accuracy for each threshold.
#' The function returns the best accuracy achieved and the corresponding purity threshold.
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
                            exclude_other_for_accuracy = FALSE,
                            other_class = "Other") {  # <--- NEW ARG

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
                                   other_class = other_class)  # <--- pass it through
  
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
      confusion_matrix <- table(updated_votes[[truth_column]], updated_votes[[out_column]])
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


#' Predict DLBCLone Classes for New Samples Using a Trained KNN Model
#'
#' Applies a previously optimized DLBCLone KNN model to predict class labels for new (test) samples.
#' This function combines the training and test feature matrices, ensures feature compatibility, and uses the
#' parameters from a DLBCLone KNN optimization run to classify the test samples. Optionally, runs in iterative mode
#' for more stable results when predicting multiple samples.
#'
#' @param train_df Data frame or matrix of features for training samples (rows = samples, columns = features).
#' @param test_df Data frame or matrix of features for test samples to be classified.
#' @param metadata Data frame with metadata for all samples, including at least a \code{sample_id} column.
#' @param core_features Optional character vector of feature names to upweight in the KNN calculation.
#' @param core_feature_multiplier Numeric. Multiplier to apply to core features (default: 1.5).
#' @param hidden_features Optional character vector of feature names to exclude from the analysis.
#' @param DLBCLone_KNN_out List. Output from a previous call to \code{DLBCLone_KNN} containing optimized parameters. (Required)
#' @param mode Character. If \code{"iterative"}, runs KNN prediction for each test sample individually (recommended for stability).
#'
#' @return A list containing the KNN prediction results for the test samples, including predicted class labels and scores.
#'
#' @details
#' - Ensures that the feature columns in \code{train_df} and \code{test_df} are compatible.
#' - If \code{mode = "iterative"}, runs KNN prediction for each test sample one at a time.
#' - Uses the parameters (e.g., k, feature weights) from the provided \code{DLBCLone_KNN_out} object.
#' - Returns the same structure as \code{DLBCLone_KNN}, with predictions for the test samples.
#'
#' @examples
#' # Assuming you have run DLBCLone_KNN to get optimized parameters:
#' # model_out <- DLBCLone_KNN(train_features, train_metadata, ...)
#' # Predict on new samples:
#' predictions <- DLBCLone_KNN_predict(
#'   train_df = train_features,
#'   test_df = new_samples,
#'   metadata = sample_metadata,
#'   DLBCLone_KNN_out = model_out
#' 
#' @export
#'
DLBCLone_KNN_predict <- function(train_df,
                                 test_df,
                                 metadata,
                                 DLBCLone_KNN_out,
                                 mode = "batch",
                                 truth_column = "lymphgen",
                                 other_class = "Other") {  # <--- NEW ARG
  if(missing(DLBCLone_KNN_out)){
    stop("DLBCLone_KNN_out must be provided, run DLBCLone_KNN first to get the optimal parameters")
  }
  nsamp = nrow(test_df)
  message(paste0("Running DLBCLone KNN individually on ", nsamp, " samples")) 
  if(nsamp > 1){
    if(mode != "iterative"){
      warning("Running DLBCLone KNN on multiple samples at once is not recommended as the result may be unstable",
              " and may not reflect the true classification of each sample. ",
              "Running in iterative mode is recommended for more stable results.")
    }
  }

  if(any(!colnames(test_df) %in% colnames(train_df))){
    stop("test_df should not contain any features that are not in train_df. ",
         "Please check the column names of test_df and train_df.")
  }
  combined_df = bind_rows(train_df, test_df)
  if(any(!colnames(train_df) %in% colnames(test_df))){
    message("filling in missing features in test_df with zeros")
    combined_df[is.na(combined_df)] = 0
  }
  if(mode == "iterative"){
    predictions_list = list()
    for(i in seq_len(nrow(test_df))){
      message("iteration:", i, "of", nrow(test_df))
      combined_df = bind_rows(train_df, test_df[i,])
      model_out = DLBCLone_KNN(
        features_df = combined_df,
        metadata = metadata,
        DLBCLone_KNN_out = DLBCLone_KNN_out,
        predict_unlabeled = TRUE,
        other_class = other_class   # <--- pass it in
      )
      predictions_list[[i]] = model_out$unlabeled_predictions
    }
    all_predictions = do.call("bind_rows", predictions_list)
    return(all_predictions)
  } else {
    model_out = DLBCLone_KNN(
      features_df = combined_df,
      metadata = metadata,
      DLBCLone_KNN_out = DLBCLone_KNN_out,
      truth_column = DLBCLone_KNN_out$truth_column,
      truth_classes = DLBCLone_KNN_out$truth_classes,
      predict_unlabeled = TRUE,
      other_class = DLBCLone_KNN_out$other_class  
    )
    return(model_out)
  }
}


#' Run DLBCLone KNN Classification
#'
#' Weighted KNN on a feature (mutation) matrix with optional upweighting of
#' user-specified "core" features, optional exclusion of "hidden" features,
#' and optional optimization of an explicit outgroup (e.g. "Other").
#'
#' This version removes hard-coded LymphGen class names and instead derives the
#' in-group classes and the outgroup column name from the arguments
#' \code{truth_classes} and \code{other_class}. It keeps backward compatibility
#' for the default LymphGen-like usage.
#'
#' @param features_df Numeric matrix/data.frame (rows = samples, cols = features).
#'                    Row names must be sample IDs.
#' @param metadata Data frame with at least \code{sample_id} and the ground-truth
#'                 label column given in \code{truth_column}.
#' @param core_features Character vector of feature names to upweight (optional).
#' @param core_feature_multiplier Numeric multiplier for \code{core_features}.
#' @param hidden_features Character vector of feature names to drop (optional).
#' @param min_k,max_k Integer K range to explore when optimizing.
#' @param truth_column Name of metadata column with ground-truth class labels.
#' @param truth_classes Character vector of all classes to consider (including
#'                      \code{other_class} if you intend to optimize for it).
#' @param other_class Name of the explicit outgroup class (default: "Other").
#' @param optimize_for_other Logical; if TRUE, computes a separate "other"
#'        score (ratio) and searches a purity threshold; if FALSE, treats all
#'        classes symmetrically.
#' @param predict_unlabeled If TRUE, re-runs KNN to classify samples that were
#'        present in \code{features_df} but not in \code{metadata}.
#' @param plot_samples Optional vector of sample_ids to keep in example plots.
#' @param DLBCLone_KNN_out Optional prior result; if supplied, its learned
#'        parameters are reused (skip optimization).
#' @param seed Random seed.
#' @param epsilon Small value added to distances before weighting.
#' @param weighted_votes If FALSE, neighbors are unweighted (equal votes).
#' @param skip_umap If TRUE, skip layout optimization plots at the end.
#'
#' @return A list with fields including:
#'   \item{predictions}{Per-sample vote/score summary and predicted labels}
#'   \item{DLBCLone_k_best_k}{Best K found}
#'   \item{DLBCLone_k_purity_threshold}{Best purity threshold (if applicable)}
#'   \item{DLBCLone_k_accuracy}{Best accuracy metric achieved}
#'   \item{truth_classes, truth_column}{Echoed arguments}
#'   \item{unlabeled_predictions}{Predictions for unlabeled samples (if requested)}
#'   \item{df}{Annotated layout for plotting (when built in this run)}
#'   \item{plot_truth, plot_predicted}{ggplots when built in this run}
#'
#' @export
DLBCLone_KNN <- function(features_df,
                         metadata,
                         core_features = NULL,
                         core_feature_multiplier = 1.5,
                         hidden_features = NULL,
                         min_k = 5,
                         max_k = 60,
                         truth_column = "lymphgen",
                         truth_classes = c("EZB", "BN2", "ST2", "MCD", "N1", "Other"),
                         other_class = "Other",
                         optimize_for_other = TRUE,
                         predict_unlabeled = FALSE,
                         plot_samples = NULL,
                         DLBCLone_KNN_out = NULL, 
                         seed = 12345, 
                         epsilon = 0.001,
                         weighted_votes = TRUE,
                         skip_umap = FALSE) {

  # In-group classes (exclude the explicit outgroup label if present)
  class_levels <- setdiff(truth_classes, other_class)

  if(!missing(DLBCLone_KNN_out)){
    core_features        <- DLBCLone_KNN_out$core_features
    core_feature_multiplier <- DLBCLone_KNN_out$core_feature_multiplier
    hidden_features      <- DLBCLone_KNN_out$hidden_features
    optimize_for_other   <- DLBCLone_KNN_out$optimize_for_other
  }

  # Upweight core features, drop hidden features
  if(!is.null(core_features)){
    if(!is.numeric(core_feature_multiplier)){
      stop("core_feature_multiplier must be a numeric value")
    }
    if(!all(core_features %in% colnames(features_df))){
      stop("core_features must be a vector of column names in features_df")
    }
    ncore <- length(core_features)
    message(paste0("multiplying ", ncore, " core features by ", core_feature_multiplier))
    features_df[, core_features] <- features_df[, core_features] * core_feature_multiplier
  }
  if(!is.null(hidden_features)){
    if(!all(hidden_features %in% colnames(features_df))){
      stop("hidden_features must be a vector of column names in features_df")
    }
    message(paste0("dropping ", length(hidden_features), " hidden features"))
    features_df <- features_df %>% dplyr::select(-dplyr::any_of(hidden_features))
  }

  # Partition rows by presence in metadata and by empty feature rows
  exclude_df  <- features_df[!row.names(features_df) %in% metadata$sample_id, , drop = FALSE]
  features_df <- features_df[ rownames(features_df) %in% metadata$sample_id, , drop = FALSE]

  df_empty <- features_df[rowSums(features_df) == 0, , drop = FALSE]
  sample_metadata_no_features <- dplyr::filter(metadata, sample_id %in% rownames(df_empty))
  features_df <- features_df[rowSums(features_df) > 0, , drop = FALSE]

  exclude_empty <- exclude_df[rowSums(exclude_df) == 0, , drop = FALSE]
  exclude_df    <- exclude_df[rowSums(exclude_df) > 0, , drop = FALSE]

  if (is.null(DLBCLone_KNN_out)) {
    # ------------------------
    # Optimize over K (and purity threshold if requested)
    # ------------------------
    overall_best_acc_k   <- 0
    overall_best_acc     <- 0
    overall_best_thresh  <- 0
    best_pred            <- NULL

    metadata <- metadata %>% dplyr::filter(sample_id %in% rownames(features_df))
    metadata_simple <- metadata %>% dplyr::select(sample_id, !!rlang::sym(truth_column))

    message("Finding all nearest neighbors up to k=", max_k + 1, " using cosine distance")
    nn_u <- uwot::umap(features_df,
                       n_neighbors = max_k + 1,
                       ret_nn = TRUE,
                       metric = "cosine",
                       seed = seed,
                       n_threads = 1,
                       batch = TRUE,
                       n_sgd_threads = 1,
                       rng_type = "deterministic")

    fkn_ids   <- nn_u$nn$cosine$idx
    fkn_dists <- nn_u$nn$cosine$dist
    rownames(fkn_dists) <- rownames(features_df)
    rownames(fkn_ids)   <- rownames(features_df)

    fkn_dists <- as.data.frame(fkn_dists) %>% dplyr::select(-1)
    if (weighted_votes) {
      fkn_weighted <- round(1 / (fkn_dists + epsilon), 7)
    } else {
      fkn_weighted <- fkn_dists
      fkn_weighted[] <- 1
    }

    fkn_ids <- as.data.frame(fkn_ids) %>% dplyr::select(-1)
    truth_index <- metadata[[truth_column]]
    names(truth_index) <- metadata$sample_id

    fkn_ids_named <- apply(fkn_ids, 2, function(x) rownames(fkn_ids)[x])
    rownames(fkn_ids_named) <- rownames(fkn_ids)

    fkn_ids_truth <- apply(fkn_ids_named, 2, function(x) truth_index[x])
    rownames(fkn_ids_truth) <- rownames(fkn_ids)
    fkn_ids_truth <- as.data.frame(fkn_ids_truth)

    fkn_ids_long <- tibble::rownames_to_column(fkn_ids_truth, "sample_id") %>%
      tidyr::pivot_longer(-sample_id, names_to = "column", values_to = "class")
    fkn_weighted_long <- tibble::rownames_to_column(fkn_weighted, "sample_id") %>%
      tidyr::pivot_longer(-sample_id, names_to = "column", values_to = "vote")

    fkn_votes <- dplyr::left_join(fkn_ids_long, fkn_weighted_long,
                                  by = c("sample_id", "column"))

    for (k in seq(max_k, min_k, by = -5)) {
      message(paste0("Running DLBCLone KNN with k=", k))

      fkn_weighted <- fkn_votes %>%
        dplyr::group_by(sample_id) %>%
        dplyr::slice_head(n = k) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(sample_id, class) %>%
        dplyr::summarize(weighted_vote = sum(vote), .groups = "drop")

      if (optimize_for_other) {
        # Treat outgroup separately; compute top group among in-groups only
        fkn_weighted <- fkn_weighted %>%
          tidyr::pivot_wider(
            id_cols     = "sample_id",
            names_from  = "class",
            values_from = "weighted_vote",
            values_fill = 0
          ) %>%
          dplyr::select(sample_id, dplyr::all_of(union(class_levels, other_class))) %>%
          dplyr::rowwise() %>%
          dplyr::mutate(
            idx             = which.max(dplyr::c_across(dplyr::all_of(class_levels))),
            top_class       = class_levels[idx],
            top_class_count = dplyr::c_across(dplyr::all_of(class_levels))[idx]
          ) %>%
          dplyr::select(-idx) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(
            top_class       = ifelse(top_class_count == 0, other_class, top_class),
            "{other_class}" := round(.data[[other_class]], 4),
            top_class_count = round(top_class_count, 4)
          ) %>%
          dplyr::rename(
            "{other_class}_score" := dplyr::all_of(other_class)
          ) %>%
          dplyr::mutate(
            by_score        = top_class,
            top_group_score = top_class_count,
            score_ratio     = round(top_group_score / .data[[paste0(other_class, "_score")]], 4)
          )

      } else {
        # All classes symmetric; still provide an {other_class}_score column
        fkn_weighted <- fkn_weighted %>%
          tidyr::pivot_wider(
            id_cols     = "sample_id",
            names_from  = "class",
            values_from = "weighted_vote",
            values_fill = 0
          ) %>%
          dplyr::select(sample_id, dplyr::all_of(union(class_levels, other_class))) %>%
          dplyr::rowwise() %>%
          dplyr::mutate(
            idx             = which.max(dplyr::c_across(dplyr::all_of(class_levels))),
            top_class       = class_levels[idx],
            top_class_count = dplyr::c_across(dplyr::all_of(class_levels))[idx]
          ) %>%
          dplyr::select(-idx) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(
            top_class = ifelse(top_class_count == 0, other_class, top_class),
            "{other_class}_score" := dplyr::if_else(
              !is.na(.data[[other_class]]), .data[[other_class]], 0
            )
          ) %>%
          dplyr::mutate(
            by_score        = top_class,
            top_group_score = top_class_count,
            score_ratio     = 10
          )
      }

      # Attach truth labels for reporting
      fkn_weighted <- dplyr::left_join(
        fkn_weighted,
        dplyr::select(metadata, sample_id, !!rlang::sym(truth_column)),
        by = "sample_id"
      )

      score_thresh <- 1.5 * k

      # Grid search purity threshold if outgroup optimization is on
      if (optimize_for_other) {
        best <- 0
        best_thresh <- 0
        for (purity_threshold in seq(3, 0, -0.05)) {
          updated_votes <- fkn_weighted %>%
            dplyr::mutate(min_score = score_thresh) %>%
            dplyr::mutate(
              by_score_opt = ifelse(
                score_ratio >= purity_threshold | top_group_score > score_thresh,
                by_score, other_class
              )
            )

          acc <- report_accuracy(updated_votes, pred = "by_score_opt", truth = truth_column)
          if (acc$mean_balanced_accuracy > best) {
            best <- acc$mean_balanced_accuracy
            best_thresh <- purity_threshold
          }
        }
      } else {
        updated_votes <- dplyr::mutate(fkn_weighted, by_score_opt = by_score)
        acc          <- report_accuracy(updated_votes, pred = "by_score_opt", truth = truth_column)
        best         <- acc$mean_balanced_accuracy
        best_thresh  <- 0
      }

      classes_for_span <- setdiff(intersect(class_levels, names(fkn_weighted)), other_class)
      if (best > overall_best_acc) {
        overall_best_acc    <- best
        overall_best_thresh <- best_thresh
        overall_best_acc_k  <- k

        updated_votes <- fkn_weighted %>%
          dplyr::mutate(min_score = score_thresh) %>%
          dplyr::mutate(DLBCLone_k  = by_score) %>%
          dplyr::mutate(DLBCLone_ko = ifelse(
            score_ratio >= best_thresh | top_group_score > score_thresh,
            by_score, other_class
          )) %>%
          dplyr::rowwise() %>%
          dplyr::mutate(
            valid_classes = {
              scores <- dplyr::c_across(dplyr::all_of(classes_for_span))
              keep_idx <- which(replace(scores > score_thresh, is.na(scores), FALSE))
              if (length(keep_idx) == 0) other_class else paste(classes_for_span[keep_idx], collapse = "/")
            }
          ) %>%
          dplyr::ungroup()

        best_pred <- updated_votes
        message(paste0(
          "New best accuracy: ", round(best, 3),
          " at k=", k, " with purity threshold: ", round(best_thresh, 2)
        ))
      } else {
        message(paste(
          "best accuracy did not improve at k=", k,
          " with purity threshold:", round(best_thresh, 2),
          " accuracy:", round(best, 3)
        ))
      }
    } # end for K

  } else {
    message("Using DLBCLone_KNN_out provided, skipping KNN run")
    best_pred            <- DLBCLone_KNN_out$predictions
    optimized_layout     <- DLBCLone_KNN_out$df
    overall_best_acc_k   <- DLBCLone_KNN_out$DLBCLone_k_best_k
    overall_best_acc     <- DLBCLone_KNN_out$DLBCLone_k_accuracy
    overall_best_thresh  <- DLBCLone_KNN_out$DLBCLone_k_purity_threshold
  }

  # ------------------------
  # Predict for unlabeled samples (optional)
  # ------------------------
  unlabeled_predictions <- data.frame()
  samples_no_metadata <- c(
    rownames(exclude_df)[!rownames(exclude_df) %in% metadata$sample_id],
    rownames(exclude_empty)[!rownames(exclude_empty) %in% metadata$sample_id]
  )

  if (predict_unlabeled && length(samples_no_metadata) > 0) {
    message("Re-running KNN to include unlabeled samples. Will use K value and thresholds from optimized model, if provided")
    message("will use newly provided features rather than recycling!")

    if (!is.null(DLBCLone_KNN_out)) {
      k <- DLBCLone_KNN_out$DLBCLone_k_best_k
      overall_best_thresh <- DLBCLone_KNN_out$DLBCLone_k_purity_threshold
      best_thresh <- overall_best_thresh
    }

    # Build an augmented metadata where the unlabeled have NA in the truth column
    placeholder_metadata <- data.frame(sample_id = samples_no_metadata, stringsAsFactors = FALSE)
    placeholder_metadata[[truth_column]] <- NA
    metadata_merge <- dplyr::bind_rows(metadata, placeholder_metadata)

    k_buffer   <- min(length(samples_no_metadata), overall_best_acc_k)
    generous_k <- overall_best_acc_k + k_buffer
    message("Finding all nearest neighbors up to k=", generous_k + 1, " using cosine distance")

    if (nrow(exclude_df) > 0) {
      features_df_merge <- dplyr::bind_rows(features_df, exclude_df)

      nn_u <- uwot::umap(features_df_merge,
                         n_neighbors = generous_k + 1,
                         ret_nn = TRUE,
                         metric = "cosine",
                         seed = seed,
                         n_threads = 1,
                         batch = TRUE,
                         n_sgd_threads = 1,
                         rng_type = "deterministic")
      fkn_ids   <- nn_u$nn$cosine$idx
      fkn_dists <- nn_u$nn$cosine$dist
      rownames(fkn_dists) <- rownames(features_df_merge)
      rownames(fkn_ids)   <- rownames(features_df_merge)

      fkn_dists <- as.data.frame(fkn_dists) %>% dplyr::select(-1)
      if (weighted_votes) {
        fkn_weighted <- round(1 / (fkn_dists + epsilon), 7)
      } else {
        fkn_weighted <- fkn_dists
        fkn_weighted[] <- 1
      }

      fkn_ids <- as.data.frame(fkn_ids) %>% dplyr::select(-1)

      truth_index <- metadata_merge[[truth_column]]
      names(truth_index) <- metadata_merge$sample_id

      fkn_ids_named <- apply(fkn_ids, 2, function(x) rownames(fkn_ids)[x])
      rownames(fkn_ids_named) <- rownames(fkn_ids)

      fkn_ids_truth <- apply(fkn_ids_named, 2, function(x) truth_index[x])
      rownames(fkn_ids_truth) <- rownames(fkn_ids)
      fkn_ids_truth <- as.data.frame(fkn_ids_truth)

      fk_neighbors_long <- as.data.frame(fkn_ids_named) %>%
        tibble::rownames_to_column("sample_id") %>%
        tidyr::pivot_longer(-sample_id, names_to = "column", values_to = "neighbor")

      fkn_ids_long <- tibble::rownames_to_column(fkn_ids_truth, "sample_id") %>%
        tidyr::pivot_longer(-sample_id, names_to = "column", values_to = "class")

      fkn_weighted_long <- tibble::rownames_to_column(fkn_weighted, "sample_id") %>%
        tidyr::pivot_longer(-sample_id, names_to = "column", values_to = "vote")

      fkn_votes <- dplyr::left_join(fkn_ids_long, fkn_weighted_long, by = c("sample_id", "column")) %>%
        dplyr::left_join(fk_neighbors_long, by = c("sample_id", "column")) %>%
        dplyr::filter(!is.na(class))

      k <- overall_best_acc_k

      fkn_weighted <- fkn_votes %>%
        dplyr::group_by(sample_id) %>%
        dplyr::slice_head(n = k) %>%
        dplyr::mutate(neighbor_id = paste(neighbor, collapse = ",")) %>%
        dplyr::ungroup()

      fkn_weighted_neighbors <- fkn_weighted %>%
        dplyr::select(sample_id, neighbor_id) %>%
        dplyr::distinct() %>%
        dplyr::filter(sample_id %in% rownames(exclude_df))

      fkn_weighted <- fkn_weighted %>%
        dplyr::group_by(sample_id, class) %>%
        dplyr::summarize(weighted_vote = sum(vote), .groups = "drop")

      if (optimize_for_other) {
        fkn_weighted <- fkn_weighted %>%
          tidyr::pivot_wider(
            id_cols     = "sample_id",
            names_from  = "class",
            values_from = "weighted_vote",
            values_fill = 0
          ) %>%
          dplyr::select(sample_id, dplyr::all_of(union(class_levels, other_class))) %>%
          dplyr::rowwise() %>%
          dplyr::mutate(
            idx             = which.max(dplyr::c_across(dplyr::all_of(class_levels))),
            top_class       = class_levels[idx],
            top_class_count = dplyr::c_across(dplyr::all_of(class_levels))[idx]
          ) %>%
          dplyr::select(-idx) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(
            top_class       = ifelse(top_class_count == 0, other_class, top_class),
            "{other_class}" := round(.data[[other_class]], 4),
            top_class_count = round(top_class_count, 4)
          ) %>%
          dplyr::rename(
            "{other_class}_score" := dplyr::all_of(other_class)
          ) %>%
          dplyr::mutate(
            by_score        = top_class,
            top_group_score = top_class_count,
            score_ratio     = round(top_group_score / .data[[paste0(other_class, "_score")]], 4)
          )
      } else {
        fkn_weighted <- fkn_weighted %>%
          tidyr::pivot_wider(
            id_cols     = "sample_id",
            names_from  = "class",
            values_from = "weighted_vote",
            values_fill = 0
          ) %>%
          dplyr::select(sample_id, dplyr::all_of(union(class_levels, other_class))) %>%
          dplyr::rowwise() %>%
          dplyr::mutate(
            idx             = which.max(dplyr::c_across(dplyr::all_of(class_levels))),
            top_class       = class_levels[idx],
            top_class_count = dplyr::c_across(dplyr::all_of(class_levels))[idx]
          ) %>%
          dplyr::select(-idx) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(
            top_class = ifelse(top_class_count == 0, other_class, top_class),
            "{other_class}_score" := dplyr::if_else(
              !is.na(.data[[other_class]]), .data[[other_class]], 0
            )
          ) %>%
          dplyr::mutate(
            by_score        = top_class,
            top_group_score = top_class_count,
            score_ratio     = 10
          )
      }

      fkn_weighted <- fkn_weighted %>%
        dplyr::filter(sample_id %in% rownames(exclude_df))

      fkn_weighted <- dplyr::left_join(
        fkn_weighted,
        dplyr::select(metadata_merge, sample_id, !!rlang::sym(truth_column)),
        by = "sample_id"
      )

      score_thresh <- 1.5 * k

      if (optimize_for_other) {
        unlabeled_predictions <- fkn_weighted %>%
          dplyr::mutate(DLBCLone_k = by_score) %>%
          dplyr::mutate(DLBCLone_ko = ifelse(
            score_ratio >= overall_best_thresh | top_group_score > score_thresh,
            by_score, other_class
          )) %>%
          dplyr::left_join(., fkn_weighted_neighbors, by = "sample_id")
      } else {
        classes_for_span <- intersect(class_levels, names(fkn_weighted))
        other_score_col  <- paste0(other_class, "_score")

        unlabeled_predictions <- fkn_weighted %>%
          dplyr::mutate(DLBCLone_k = by_score) %>%
          dplyr::mutate(DLBCLone_ko = by_score) %>%
          dplyr::rowwise() %>%
          dplyr::mutate(
            valid_classes = {
              scores   <- c_across(all_of(classes_for_span))
              keep_idx <- which(replace(scores > .data[[other_score_col]], is.na(scores), FALSE))
              if (length(keep_idx) == 0) other_class else paste(classes_for_span[keep_idx], collapse = "/")
            }
          ) %>%
          dplyr::ungroup() %>%
          dplyr::left_join(., fkn_weighted_neighbors, by = "sample_id")
      }
    } else {
      other_score_col <- paste0(other_class, "_score")
      unlabeled_predictions <- data.frame(sample_id = rownames(exclude_empty), stringsAsFactors = FALSE)
      unlabeled_predictions[[truth_column]] <- NA
      unlabeled_predictions <- unlabeled_predictions %>%
        dplyr::mutate(
          DLBCLone_k   = other_class,
          DLBCLone_ko  = other_class,
          by_score     = NA,
          top_group_score = NA,
          score_ratio  = NA
        )
      # Ensure truth_classes columns exist as NA
      for (cls in truth_classes) {
        if (!cls %in% colnames(unlabeled_predictions)) {
          unlabeled_predictions[[cls]] <- NA
        }
      }
      unlabeled_predictions[[other_score_col]] <- NA
      # no neighbor_id information in this branch
    }
  }

  # ------------------------
  # Assemble return object
  # ------------------------
  if (is.null(DLBCLone_KNN_out)) {

    if (nrow(sample_metadata_no_features) > 0) {
      sample_metadata_no_features <- sample_metadata_no_features %>%
        dplyr::select(sample_id, !!rlang::sym(truth_column)) %>%
        dplyr::mutate(
          DLBCLone_k = NA,
          DLBCLone_ko = other_class,
          by_score = NA,
          top_group_score = NA,
          score_ratio = NA
        )
      best_pred <- dplyr::bind_rows(best_pred, sample_metadata_no_features)
    }

    format_for_output <- function(x) {
      x_rounded <- round(x, 4)
      ifelse(
        is.na(x_rounded),
        NA_character_,
        ifelse(
          x_rounded == 0,
          "0",
          sub("\\.?0+$", "", format(x_rounded, scientific = FALSE, trim = TRUE, nsmall = 0))
        )
      )
    }

    to_return <- list(
      predictions = best_pred %>%
        dplyr::mutate(dplyr::across(where(is.numeric), ~ format_for_output(.))),
      DLBCLone_k_best_k        = overall_best_acc_k,
      DLBCLone_k_purity_threshold = overall_best_thresh,
      DLBCLone_k_accuracy      = overall_best_acc,
      truth_classes            = truth_classes,
      truth_column             = truth_column,
      sample_metadata_no_features = sample_metadata_no_features,
      core_feature_multiplier  = core_feature_multiplier,
      core_features            = core_features,
      hidden_features          = hidden_features,
      seed                     = seed,
      optimize_for_other       = optimize_for_other,
      unlabeled_predictions    = NULL
    )
  } else {
    to_return <- DLBCLone_KNN_out
  }

  # Fill in neighbors / predictions for unlabeled (if requested)
  to_return$unlabeled_neighbors   <- NULL
  if (predict_unlabeled && length(samples_no_metadata) > 0) {
    if (exists("fkn_weighted_neighbors")) {
      to_return$unlabeled_neighbors <- fkn_weighted_neighbors %>%
        tidyr::separate(neighbor_id, into = paste0("N", c(1:k)), sep = ",")
    }
    to_return$unlabeled_predictions <- unlabeled_predictions

    if (is.null(DLBCLone_KNN_out) & !skip_umap) {
      message("Optimizing graph layout for visualization, predict_unlabeled = TRUE")
      df_show <- dplyr::bind_rows(features_df, exclude_df)

      optimized <- make_and_annotate_umap(df_show, metadata = metadata)$df
      optimized <- optimized %>%
        dplyr::left_join(
          dplyr::bind_rows(
            dplyr::select(best_pred, -!!rlang::sym(truth_column)),
            dplyr::select(unlabeled_predictions, -!!rlang::sym(truth_column))
          ),
          by = "sample_id"
        ) %>%
        dplyr::mutate(!!rlang::sym(truth_column) := ifelse(is.na(.data[[truth_column]]), DLBCLone_ko, .data[[truth_column]]))
    } else {
      optimized <- if (!is.null(DLBCLone_KNN_out)) DLBCLone_KNN_out$features_df else NULL
    }
  } else {
    if (is.null(DLBCLone_KNN_out)) {
      message("Optimizing graph layout for visualization")
      optimized <- make_and_annotate_umap(features_df, metadata = metadata)$df
      optimized <- optimized %>%
        dplyr::left_join(dplyr::select(best_pred, -!!rlang::sym(truth_column)), by = "sample_id") %>%
        dplyr::mutate(!!rlang::sym(truth_column) := ifelse(is.na(.data[[truth_column]]), DLBCLone_ko, .data[[truth_column]]))
    } else {
      message("Using DLBCLone_KNN_out provided, skipping graph layout optimization")
      optimized <- DLBCLone_KNN_out$features_df
    }
  }

  if(!is.null(optimized) & is.null(DLBCLone_KNN_out)){
    to_return$plot_truth     <- basic_umap_scatterplot(optimized, plot_samples, colour_by = truth_column)
    to_return$plot_predicted <- basic_umap_scatterplot(optimized, plot_samples, colour_by = "DLBCLone_ko")
    to_return$df             <- optimized
  }

  # Ensure features_df is included in the return for downstream plotting
  if (predict_unlabeled && length(samples_no_metadata) > 0) {
    to_return$features_df <- if (exists("features_df_merge")) features_df_merge else features_df
  } else {
    to_return$features_df <- features_df
  }

  to_return$type <- "DLBCLone_KNN"
  to_return$pred_column <- "DLBCLone_ko"
  to_return$other_class = other_class
  return(to_return)
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
#' @import caret
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
                           truth_column = "lymphgen",
                           optimize_for_other = FALSE,
                           
                           eval_group = NULL,
                           min_k=3,
                           max_k=23,
                           verbose = FALSE,
                           seed = 12345,
                           maximize = "balanced_accuracy", #or "harmonic_mean" or "accuracy"
                           exclude_other_for_accuracy = FALSE,
                           weights_opt = c(TRUE)
                           ) {
  macro_f1 <- function(truth, pred, drop_other = TRUE, other_label = "Other", na_rm = TRUE) {
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

  exclude_other_for_accuracy = TRUE

  na_opt = c("drop")
  num_class = length(truth_classes)
  
  threshs = seq(0,0.9,0.1)

  
  ks = seq(min_k,max_k,2)
  results <- data.frame()
  best_params <- data.frame()

  use_w = TRUE
  best_acc = 0
  best_acc_w = 0
  best_pred_w = NULL
  best_k_w = 0
  best_fit = NULL
  this_accuracy = 0
  best_pred = NULL
  other_pred = NULL
  best_w_score_thresh = NULL
  best_w_purity = NULL
  w_best_pred = NULL
  ignore_self = TRUE

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
      if(!missing(combined_mutation_status_df)){
        message("ignoring mutation status data frame and using features from umap_out instead")
      }
      outs = make_and_annotate_umap(df=umap_out$features,
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
    

    for(use_w in weights_opt){
      for(k in ks){
        best_w_acc = NULL
        message(paste("K:",k))
        test_coords = filter(outs$df,!!sym(truth_column) %in% truth_classes) %>% select(V1,V2)
        train_coords = filter(outs$df,!!sym(truth_column) %in% truth_classes) %>% select(V1,V2)
        train_labels = filter(outs$df,!!sym(truth_column) %in% truth_classes) %>% pull(!!sym(truth_column))
        train_ids = filter(outs$df,!!sym(truth_column) %in% truth_classes) %>% pull(sample_id)
        test_ids = filter(outs$df,!!sym(truth_column) %in% truth_classes) %>% pull(sample_id)
        rownames(train_coords) = train_ids
        rownames(test_coords) = test_ids



        pred = weighted_knn_predict_with_conf(
          train_coords = train_coords,
          train_labels = train_labels,
          test_coords = test_coords,
          k=k,
          conf_threshold =threshold,
          other_class="Other",
          use_weights = use_w,
          ignore_self = ignore_self,
          verbose = verbose)

        for(threshold in threshs){
          if(verbose){
            print(paste("k:",k,"threshold:",threshold,"use_weights:",use_w,"na_option:",na_option))
          }
          pred_thresh = mutate(pred,
            predicted_label = ifelse(confidence >= threshold, predicted_label, "Other"))
          
          if(is.null(best_w_acc)){
            #DLBCLone_wo
            some_outs = filter(outs$df,!!sym(truth_column) %in% truth_classes)
            pred_with_truth = bind_cols(some_outs ,pred)
            if(verbose){
              print(dim(pred_with_truth))
              print(table(pred_with_truth[[truth_column]]))
            }

            pred_w = 
             optimize_purity(
              vote_df = pred_with_truth,
              mode = "DLBCLone_w",
              truth_column = truth_column,
              all_classes = truth_classes,
              optimize_by = maximize,
              k = k
             )
             
            if(verbose){
              print("running optimize_purity for the first time at k:", k)
            }
          
            best_w_acc = pred_w$best_accuracy
            if(best_w_acc > best_acc_w){
              
              best_acc_w = best_w_acc
              best_k_w = k
              best_fit = pred_with_truth
              best_pred_w = pred_w$predictions
              best_w_purity = pred_w$best_purity_threshold
              best_w_score_thresh = pred_w$score_thresh
              #print(head(pred_w$predictions))
              reported_accuracy = report_accuracy(pred_w$predictions,pred = "DLBCLone_wo")
              conc = reported_accuracy$accuracy_no_other
              f1 = macro_f1(pred_w$predictions[[truth_column]],
                            pred_w$predictions$DLBCLone_wo,
                            drop_other = T)
              new_f1 = reported_accuracy$macro_f1
              harmonic_mean = reported_accuracy$harmonic_mean
              mba = reported_accuracy$mean_balanced_accuracy
              cr = reported_accuracy$classification_rate
                
              message("new best DLBCLone_wo:")
              message(paste(
                "Concordance:",
                round(conc,3),
                "F1:",
                round(f1,5),
                "MBA:",
                round(best_w_acc,3),
                "Classification rate:",
                round(cr,5),
                "Harmonic mean:",
                round(harmonic_mean,5),
                "purity:",
                pred_w$best_purity_threshold,
                "Score:",
                best_w_score_thresh
              ))
            }
            best_pred_w = pred_w$predictions
            
          }


          
          train_d = filter(outs$df,!!sym(truth_column) %in% truth_classes)
          xx_d = bind_cols(train_d ,pred_thresh)

          if("Other" %in% truth_classes){
            xx_d[[truth_column]] = factor(xx_d[[truth_column]])
          }else{
            test_coords = filter(outs$df,!!sym(truth_column) %in% "Other" | is.na(!!sym(truth_column))) %>% select(V1,V2)
            train_coords = filter(outs$df,!!sym(truth_column) %in% truth_classes) %>% select(V1,V2)
            train_labels = filter(outs$df,!!sym(truth_column) %in% truth_classes) %>% pull(!!sym(truth_column))
            n_other = nrow(test_coords)

            if(!"Other" %in% truth_classes && n_other > 0){
              pred_other = weighted_knn_predict_with_conf(
                train_coords = train_coords,
                train_labels = train_labels,
                test_coords = test_coords,
                k=k,
                conf_threshold =threshold,
                other_class="Other",
                use_weights = use_w,
                ignore_self = ignore_self,
                verbose = verbose)

              xx_o = bind_cols(filter(outs$df,!!sym(truth_column) == "Other" | is.na(!!sym(truth_column))) ,pred_other)

            }
            xx_d[[truth_column]] = factor(xx_d[[truth_column]],levels = c(unique(xx_d[[truth_column]]),"Other"))

          }
          true_factor = xx_d[[truth_column]]
          pred_factor = factor(xx_d$predicted_label,levels = levels(xx_d[[truth_column]]))

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
            print(conf_matrix$overall)
            print(sn)
          }
          
          overall_accuracy <- conf_matrix$overall[["Accuracy"]]
          
          overall_sensitivity<- mean(sn[!names(sn) == "Class: Other"], na.rm = TRUE)
  
          if(optimize_for_other){
   
            optimized_accuracy_and_thresh = optimize_outgroup(pred_factor,
                                             true_factor,
                                             xx_d$other_score,
                                             all_classes = truth_classes,
                                             maximize = maximize,
                                             exclude_other_for_accuracy = exclude_other_for_accuracy)
            
            out_opt_thresh = optimized_accuracy_and_thresh$threshold
            out_opt_acc = optimized_accuracy_and_thresh$average_accuracy
            

          }else{
            out_opt_acc = 0
            out_opt_thresh = 0
          }
          if(exclude_other_for_accuracy){
            mean_balanced_accuracy = mean(bal_acc[!names(bal_acc) == "Class: Other"], na.rm = TRUE)
          }else{
            mean_balanced_accuracy = mean(bal_acc)
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
          
            xx_d = mutate(xx_d, predicted_label_optimized = ifelse(other_score > out_opt_thresh, "Other", predicted_label))
            no_other_acc = report_accuracy(xx_d,pred = "predicted_label_optimized")$mean_balanced_accuracy
            message("new best accuracy DLBCLone_io:")
            message(paste(
                          best_acc, 
                          "without other: ",
                          no_other_acc, 
                          "sensitivity:",
                          overall_sensitivity,
                          "threshold_outgroup:",
                          row$threshold_outgroup))

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
  best_params$num_features = ncol(umap_out$features)
  best_params$seed = seed

  test_coords = outs$df %>% select(V1,V2)

  train_coords = outs$df %>% select(V1,V2)
  train_ids = outs$df %>% pull(sample_id)
  rownames(train_coords) = train_ids
  train_labels = outs$df %>% pull(!!sym(truth_column))
  rownames(test_coords) = outs$df %>% pull(sample_id)
  

  if(verbose){
    print(paste("TOP score threshold:",best_w_score_thresh, "purity:", best_w_purity))  
  }
  

  pred = weighted_knn_predict_with_conf(
            train_coords = train_coords,
            train_labels = train_labels,
            test_coords = test_coords,
            k=best_k_w,# re-run with best k for weighted voting
            #conf_threshold =best_params$threshold,
            use_weights = best_params$use_weights,
            ignore_self = ignore_self,
            verbose = verbose,
            track_neighbors = TRUE)
  
  pred_with_truth_full = bind_cols(outs$df %>% select(sample_id, !!sym(truth_column)), pred)

  best_pred_w = process_votes(df=pred_with_truth_full,
                              group_labels=truth_classes,
                              k=best_k_w) 
  
  best_pred_w = best_pred_w %>%
                mutate(DLBCLone_w = ifelse(score_ratio >= best_w_purity | 
                    top_group_score > best_w_score_thresh,
                                           by_score,
                                           "Other")) 

  
  if(optimize_for_other){
    
    pred = mutate(pred,predicted_label_optimized = ifelse(other_score > best_params$threshold_outgroup,
                                                          "Other",
                                                          predicted_label))
  }else{
    pred = mutate(pred,predicted_label_optimized = predicted_label)
  }

  #new naming convention
  pred = mutate(pred, DLBCLone_i = predicted_label, DLBCLone_io = predicted_label_optimized)
  best_pred_w = rename(best_pred_w, DLBCLone_wo = DLBCLone_w) #optimized
  best_pred_w = rename(best_pred_w, DLBCLone_w = by_score) #greedy
  #check accuracy again
  print(dim(pred))
  xx_d = bind_cols(outs$df,pred)
  print(dim(xx_d))
  xx_d = left_join(xx_d, select(best_pred_w, sample_id, score_ratio, top_group_score, DLBCLone_w, DLBCLone_wo), by="sample_id")
  acc_check_w = report_accuracy(xx_d,pred = "DLBCLone_w")
  acc_check_wo = report_accuracy(xx_d,pred = "DLBCLone_wo")
  message(paste("DLBCLone_w accuracy:",
  round(acc_check_w$mean_balanced_accuracy,3), 
  "DLBCLone_wo accuracy:",
  round(acc_check_wo$mean_balanced_accuracy,3),
  "DLBCLone_w Concordance:", round(acc_check_w$no_other,3),
  "DLBCLone_wo Concordance:", round(acc_check_wo$no_other,3)))
  to_ret = list(params=results,
                best_params = best_params,
                model=umap_out$model,
                features=umap_out$features,
                k_DLBCLone_i = best_params$k,
                threshold_DLBCLone_i = best_params$threshold,
                theshold_outgroup_DLBCLone_i = best_params$threshold_outgroup,
                k_DLBCLone_w = best_k_w,
                purity_DLBCLone_w = best_w_purity,
                score_thresh_DLBCLone_w = best_w_score_thresh,
                df=outs$df, 
                predictions=xx_d)
  if(!"Other" %in% truth_classes && n_other > 0){
    to_ret[["predictions_other"]] = xx_o
    to_ret[["predictions_combined"]] = bind_rows(xx_o,best_pred)
  }
  #TODO: roll components of umap_out into to_ret
  if(any(!names(outs) %in% names(to_ret))){
    missing_names = names(outs)[!names(outs) %in% names(to_ret)]
    for(mn in missing_names){
      to_ret[[mn]] = outs[[mn]] 
    }
  }
  to_ret$best_pred_w = best_pred_w
  to_ret$truth_classes = truth_classes
  to_ret$optimize_for_other = optimize_for_other
  to_ret$truth_column = truth_column
  to_ret$type = "DLBCLone_optimize_params"
  return(to_ret)
}



#' Weighted k-nearest neighbor with confidence estimate
#'
#' @param train_coords Data frame of coordinates for labeled (training) samples.
#'        One row per sample, columns are features (typically UMAP V1, V2).
#' @param train_labels Character/factor vector of labels for training samples.
#' @param test_coords  Data frame of coordinates for samples to classify
#'        (same columns/space as train_coords).
#' @param k Integer; number of neighbors to consider.
#' @param epsilon Numeric; small value added to distances before weighting
#'        (when use_weights = TRUE). Default: 0.1.
#' @param conf_threshold Optional numeric; minimum confidence for assigning a
#'        class. If provided and confidence < threshold, sample is assigned
#'        \code{other_class}.
#' @param other_class Name of the outgroup class to treat specially when
#'        \code{separate_other = TRUE}. Default: "Other".
#' @param verbose Logical; print verbose info. Default: FALSE.
#' @param use_weights Logical; inverse-distance weights (1 / (d + epsilon)).
#'        If FALSE, neighbors contribute equally. Default: TRUE.
#' @param ignore_self Logical; drop a zero-distance self-neighbor. Default: TRUE.
#' @param track_neighbors Logical; append neighbor diagnostics to output.
#'        Default: TRUE.
#' @param separate_other Logical; when TRUE, exclude neighbors labeled
#'        \code{other_class} from the main weighted vote and report their
#'        influence separately (as \code{other_*} columns). Default: TRUE.
#' @param max_neighbors Integer; maximum neighbors to retrieve from the
#'        search index before trimming to k. Default: 500.
#'
#' @return Data frame with rows = test samples and columns:
#'   \item{predicted_label}{the predicted class}
#'   \item{confidence}{predicted class weight / total weight}
#'   If \code{track_neighbors = TRUE}, additional columns:
#'   \item{other_score}{relative weight of outgroup vs predicted class}
#'   \item{neighbor_id}{comma-separated neighbor sample IDs}
#'   \item{neighbor}{comma-separated neighbor indices (in train order)}
#'   \item{distance}{comma-separated neighbor distances}
#'   \item{label}{comma-separated neighbor labels (in-group only if separate_other=TRUE)}
#'   \item{vote_labels}{comma-separated unique labels contributing to weights}
#'   \item{weighted_votes}{comma-separated weights per \code{vote_labels}}
#'   \item{neighbors_other}{count of outgroup neighbors closer than the farthest in-group neighbor}
#'   \item{other_weighted_votes}{sum of outgroup weights closer than the farthest in-group neighbor}
#'   \item{total_w}{sum of weights for in-group neighbors used}
#'   \item{pred_w}{weight supporting the predicted class}
#'
#' @import FNN
#' @export
#' 
#' 
weighted_knn_predict_with_conf <- function(train_coords,
                                           train_labels,
                                           test_coords,
                                           k,
                                           epsilon = 0.1,
                                           conf_threshold = NULL,
                                           other_class = "Other",
                                           verbose = FALSE,
                                           use_weights = TRUE,
                                           ignore_self = TRUE,
                                           track_neighbors = TRUE,
                                           separate_other = TRUE,
                                           max_neighbors = 500) { 
  
  if (nrow(train_coords)==0 || nrow(test_coords) == 0) {
    print("train_coords:")
    print(nrow(train_coords))
    print("test:")
    print(nrow(test_coords))
    stop("train_coords and test_coords must be data frames with at least one row")
  }
  # get the 100 nearest neighbors
  nn <- get.knnx(train_coords, test_coords, max_neighbors)
  all_neighbors = data.frame()
  preds <- character(nrow(test_coords))
  confs <- numeric(nrow(test_coords))

  train_labels = as.character(train_labels)
  for (i in 1:nrow(test_coords)) {

    self = rownames(test_coords)[i]

    if(verbose){
      print(paste("index:",i))
      print(test_coords[i,])
    }
    neighbors <- nn$nn.index[i, ]
    distances <- nn$nn.dist[i, ]
    if(ignore_self){
      nns = length(neighbors) 
      neighbor_ids = rownames(train_coords)[neighbors]
      
      neighbors <- neighbors[neighbor_ids != self]
      distances <- distances[neighbor_ids != self]
      neighbor_ids = neighbor_ids[neighbor_ids != self]
      
      neighbor_labels <- train_labels[neighbors]
      
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
    other_mask  <- (neighbor_labels == other_class)
    other_dists <- distances[other_mask]
    other_ids <- neighbor_ids[other_mask]
    other_labels <- neighbor_labels[other_mask]
    valid <- !is.na(neighbor_labels)
    if(!all(other_labels == other_class)){
      stop("logic error: not all other_labels are other_class")
    }
    #num_other_neighbors = sum(neighbor_labels == "Other")
    #other_dists = distances[neighbor_labels == "Other"]
    if (separate_other) valid <- valid & !other_mask

    neighbor_labels <- neighbor_labels[valid]
    weights <- weights[valid]
    distances <- distances[valid]
    neighbors <- neighbors[valid]
    #number of neighbours should be, at least, k - 1. If less than that, warn the user
    if(length(neighbors) < k-1){
      print(paste("Warning: number of neighbors is less than k-1."))
      print(paste("i:", i,"k:",k))
      print(paste("num_neighbors:",length(neighbors)))
      print(table(valid))
    }
    #now take the first k neighbors
    if(length(neighbor_labels) > k){
      neighbor_labels = neighbor_labels[1:k]
      weights = weights[1:k]
      distances = distances[1:k]
      neighbors = neighbors[1:k]
    }

    others_closer = which(other_dists < max(distances))
    
    others_distances = other_dists[others_closer]
    other_ids = other_ids[others_closer]
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
          other_neighbor = paste(other_ids,collapse=","),
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
      preds[i] <- other_class
      confs[i] <- NA
      total_weight  <- 0
      pred_weight <- 0
    } else {
      #print(weighted_votes)
      predicted_label <- names(which.max(weighted_votes))
      total_weight <- sum(weighted_votes)
      pred_weight <- weighted_votes[predicted_label]
      confidence <- pred_weight / total_weight

      # Confidence thresholding
      #if (!is.null(conf_threshold) && confidence < conf_threshold) {
      #  preds[i] <- na_label
      #  confs[i] <- confidence
      #} else {
      preds[i] <- predicted_label
      confs[i] <- confidence
      #}
    }

    if(track_neighbors){
      # Create a data frame to store neighbors, distances, and weights
      if(separate_other){
        rel_other = other_weighted_votes / pred_weight
      }else{
        rel_other = 0
      }

      neighbor_info <- data.frame(
        other_score = rel_other,
        neighbor_id = paste(rownames(train_coords)[neighbors],collapse=","),
          neighbor = paste(neighbors,collapse=","),
          distance = paste(round(distances, 3),collapse=","),
          label = paste(neighbor_labels,collapse=","),
        other_neighbor = paste(other_ids,collapse=","),
        vote_labels = paste(names(weighted_votes),collapse=","),
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
    if(nrow(to_return ) != nrow(test_coords)){
      print("mismatch in row number for to_return and test_coords")
      print(nrow(to_return))
      print(nrow(test_coords))
      print(head(to_return))
      print(head(test_coords))
      stop("")
    }
    rownames(to_return) <- rownames(test_coords)
  }
  return(to_return)
}

#' @export
prepare_single_sample_DLBCLone <- function(optimized_model,seed=12345){
  if(!optimized_model$type == "DLBCLone_optimize_params"){
    stop("Input must be the output of predict_single_sample_DLBCLone")
  }
  if(!  "projection" %in% names(optimized_model)){
      projection <- make_and_annotate_umap(
        df = optimized_model$features,
        umap_out = optimized_model,
        ret_model = FALSE,
        seed = seed,
        join_column = "sample_id",
        na_vals = optimized_model$best_params$na_option
      )
    }else{
      stop("Model already contains a projection. Nothing to do!")
    }
  optimized_model$projection = projection
  return(optimized_model)
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
#' @param ignore_self Set to TRUE to avoid considering a neighbor with the same ID
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
#' @returns a list of data frames with the predictions, the UMAP input, the model, and a ggplot object
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
    projection,
    umap_out,
    best_params,
    optimized_model = NULL,
    other_df,
    ignore_self = FALSE,
    truth_classes = c("EZB","MCD","ST2","N1","BN2","Other"),
    drop_unlabeled_from_training=TRUE,
    make_plot = TRUE,
    annotate_accuracy = FALSE,
    label_offset = 2,
    title1="GAMBL",
    title2="predicted_class_for_HighConf",
    title3 ="predicted_class_for_Other",
    seed = 12345,
    max_neighbors = 500,
    other_class = "Other"
){
    set.seed(seed)
    if(is.null(optimized_model)){
      warning("optimized_model will become a required argument in the future. Please update your code accordingly")
    }
    if(ignore_self){
        # Allow overlapping samples: rename test duplicates temporarily
        dupes <- intersect(train_df$sample_id, test_df$sample_id)
        if(length(dupes) > 0){
            test_df <- test_df %>%
                mutate(sample_id = ifelse(
                    sample_id %in% dupes,
                    paste0(sample_id, "_test"),
                    sample_id
                ))
        }
    } else {
        # Drop overlaps to prevent rowname collisions
        dupes <- intersect(train_df$sample_id, test_df$sample_id)

    }


    train_df = train_df %>%
        column_to_rownames("sample_id") %>%

        rownames_to_column("sample_id")
    train_id <- train_df$sample_id

    test_df = test_df %>%
        column_to_rownames("sample_id") %>%

        rownames_to_column("sample_id")
    test_id <- test_df$sample_id
    
    combined_df <- bind_rows(train_df, test_df)
    
    if(missing(projection)){
      projection <- make_and_annotate_umap(
        df = combined_df,
        umap_out = umap_out,
        ret_model = FALSE,
        seed = seed,
        join_column = "sample_id",
        na_vals = best_params$na_option
      )
    }else{
      message("Using provided projection for prediction")
    }
    

    train_coords = dplyr::filter(
        projection$df,
        sample_id %in% train_id
    ) %>% 
        select(sample_id,V1,V2) %>%
        column_to_rownames("sample_id")
    print("TRAIN:")
    print(dim(train_coords))
    train_df_proj = dplyr::filter(
        projection$df,
        sample_id %in% train_id
    ) %>% 
    select(sample_id,V1,V2) %>%
        left_join( #Join to the incoming metadata rather than trusting the metadata in the projection
            train_metadata %>% select(
                sample_id, 
                lymphgen
            ), 
            by = "sample_id"
        )

    train_labels = train_df_proj %>%
        pull(lymphgen) 

    #Obtain UMAP coordinates for the test sample(s)
    print(paste("num rows:",nrow(test_df)))
    test_projection <- make_and_annotate_umap(
        df = test_df,
        umap_out = umap_out,
        ret_model = FALSE,
        seed = seed,
        join_column = "sample_id",
        na_vals = best_params$na_option
      )
    test_coords = dplyr::filter(
        test_projection$df,
        sample_id %in% test_id
    ) %>% 
        select(sample_id,V1,V2) %>%
        column_to_rownames("sample_id")

    predict_training = FALSE
    if(predict_training){
      train_pred = weighted_knn_predict_with_conf(
        train_coords = train_coords,
        train_labels = train_labels,
        test_coords = train_coords, # <- predicitng training on self
        k = best_params$k,
        conf_threshold = best_params$threshold,
        other_class = other_class,
        use_weights = best_params$use_weights,
        ignore_self = ignore_self
      )

      train_pred = rownames_to_column(train_pred, var = "sample_id")
    }
    
 
    test_pred = weighted_knn_predict_with_conf(
        train_coords = train_coords,
        train_labels = train_labels,
        test_coords = test_coords,
        k = best_params$k,
        conf_threshold = best_params$threshold,
        other_class = other_class,
        use_weights = best_params$use_weights,
        ignore_self = ignore_self,
        max_neighbors = max_neighbors
    )

    test_pred = rownames_to_column(test_pred, var = "sample_id")
    #test_pred = bind_cols(test_pred, test_coords)
    if(!is.null(optimized_model)){
      test_pred = mutate(test_pred, 
                       DLBCLone_i = predicted_label,
                       DLBCLone_io = ifelse(other_score > best_params$threshold_outgroup,
                                                          other_class,
                                                          predicted_label)) 

    }else{
      test_pred = mutate(test_pred, 
                       DLBCLone_i = predicted_label,
                       DLBCLone_io = ifelse(other_score > best_params$threshold_outgroup,
                                                         other_class,
                                                          predicted_label)) 
    }
  
    anno_umap = select(test_projection$df, sample_id, V1, V2) 

    anno_out = left_join(test_pred,anno_umap,by="sample_id") %>%
        mutate(label = paste(sample_id,predicted_label,round(confidence,3)))

    anno_out = anno_out %>%
    mutate(
        V1 = as.numeric(V1),
        V2 = as.numeric(V2),
        label = as.character(label)
    )
    if(predict_training){
      predictions_train_df = left_join(train_pred, projection$df, by = "sample_id") 
    }else{
      predictions_train_df = filter(projection$df, sample_id %in% train_id) %>%
        select(sample_id, V1, V2) 
    }
    #This had assumed that projection$df contains the test samples!
    #predictions_test_df = left_join(test_pred, projection$df, by = "sample_id")
    predictions_test_df = left_join(test_pred, test_projection$df, by = "sample_id") 

    predictions_df = bind_rows(predictions_train_df %>% select(sample_id, V1, V2), 
                               predictions_test_df  %>% select(sample_id, V1, V2))

    if(!is.null(optimized_model)){
      best_w_purity = optimized_model$best_w_purity
      best_w_score_thresh = optimized_model$best_w_score_thresh
      
      predictions_test_df = process_votes(
        df=predictions_test_df,
        group_labels=optimized_model$truth_classes,
        k=optimized_model$k_DLBCLone_w) 
      
      predictions_test_df = predictions_test_df %>%
        mutate(DLBCLone_w = by_score,
          DLBCLone_wo = ifelse(score_ratio >= optimized_model$purity_DLBCLone_w | 
                                  top_group_score > optimized_model$score_thresh_DLBCLone_w, 
                                by_score, 
                                other_class))

    }
    to_return = list(
        prediction = predictions_test_df, 
        projection = projection$df,
        umap_input = umap_out$features, 
        model=umap_out$model,
        features_df = combined_df %>% column_to_rownames("sample_id"),
        df = predictions_df,
        anno_df = predictions_df %>% left_join(.,train_metadata,by="sample_id"),
        anno_out = anno_out,
        type = "predict_single_sample_DLBCLone"
      )

    return(to_return)
}

#' model storage for DLBCLone outputs
#' 
#' @param optimized_params List containing the optimized parameters from DLBCLone_optimize_params
#' @param path Path to save the files
#' @param name_prefix Prefix for the saved files, all files will be in path and start with name_prefix
#'
#' @returns saves the files to the specified path
#' 
#' @import uwot
#' @import readr
#' 
#' @export
#'
#' @examples
#' DLBCLone_save_optimized(
#'  optimzied_params = optimized_params,
#'  path="/save_optimized/trial_folder",
#'  name_prefix="test_A"
#' )
#'

DLBCLone_save_optimized = function( 
    optimized_params=NULL,
    path="models/",
    name_prefix="test"
){
  # all files will be in path and start with name_prefix
  prefix = paste0(path,"/",name_prefix)
  
  out_mut = paste0(prefix,"_mutation_status_df.tsv")
  write_tsv(optimized_params$features,file=out_mut)

  out_meta = paste0(prefix,"_metadata.tsv")
  write_tsv(optimized_params$df,file=out_meta)
  
  out_param = paste0(prefix,"_optimized_best_params.rds")
  saveRDS(optimized_params$best_params,file=out_param)
  
  out_model = paste0(prefix,"_optimized_uwot.rds")
  save_uwot(optimized_params$model,file=out_model)

  out_pred = paste0(prefix,"_optimized_pred.tsv")
  write_tsv(optimized_params$predictions,file=out_pred)

  out_classes = paste0(prefix,"_classes.txt")
  write.table(optimized_params$truth_classes,file=out_classes,quote=F,row.names=F)
}


#' load previously saved DLBCLone model and parameters
#' 
#' @param path Path to open saved the files
#' @param name_prefix Prefix of the saved files, all files will be in path and start with name_prefix
#'
#' @returns saves the files to the specified path
#' 
#' @import uwot
#' @import readr
#' @import dplyr
#' 
#' @export
#'
#' @examples
#' load_optimized <- DLBCLone_load_optimized(
#'   path="/save_optimized/trial_folder",
#'   name_prefix="test_A"
#' )
#'

DLBCLone_load_optimized <- function( # set sample_id to rownames
  path="models/",
  name_prefix="test"
){
  #all files will be in path and start with name_prefix
  prefix = paste0(path,"/",name_prefix)
  
  load_mut = paste0(prefix,"_mutation_status_df.tsv")
  load_meta = paste0(prefix,"_metadata.tsv")
  load_classes = paste0(prefix,"_classes.txt")
  load_param = paste0(prefix,"_optimized_best_params.rds")
  load_model = paste0(prefix,"_optimized_uwot.rds")
  load_pred = paste0(prefix,"_optimized_pred.tsv")

  required_files <- c(load_mut, load_meta, load_classes, load_param, load_model, load_pred)

  # Check existence
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0) {
    stop("The following required files are missing:\n", paste(missing_files, collapse = "\n"))
  }

  mut_df <- read_tsv(load_mut)
  metadata <- read_tsv(load_meta) %>%
    mutate(lymphgen = as.factor(lymphgen))
  classes <- read.table(load_classes)
  best_params <- readRDS(load_param)
  uwot_model <- load_uwot(load_model)
  predictions <- read_tsv(load_pred)

  return(list(
    df = mut_df,
    metadata = metadata,
    truth_classes = classes,
    best_params = best_params,
    model = uwot_model,
    predictions = predictions
  ))
}

#' Train a Gaussian Mixture Model for DLBCLone Classification
#'
#' Fits a supervised Gaussian mixture model (GMM) to UMAP-projected data for DLBCLone subtypes, excluding samples labeled "Other".
#' Assigns class predictions and optionally reclassifies samples as "Other" based on probability and density thresholds.
#'
#' @param umap_out List. Output from \code{make_and_annotate_umap}, containing a data frame with UMAP coordinates and truth labels.
#' @param probability_threshold Numeric. Minimum posterior probability required to assign a class (default: 0.5).
#' @param density_max_threshold Numeric. Minimum maximum density required to assign a class (default: 0.05).
#' @param cohort Optional character. Cohort label to annotate predictions.
#'
#' @details
#' - Uses \code{MclustDA} to fit a supervised mixture model to the UMAP coordinates (V1, V2) and class labels.
#' - Predicts class membership and computes per-class densities for each sample.
#' - Samples with low maximum probability or density are reclassified as "Other".
#' - Returns both raw and thresholded class assignments, respectively under the columns DLBCLone_g and DLBCLone_go.
#'
#' @return A list with:
#'   \item{gaussian_mixture_model}{Fitted \code{MclustDA} model object}
#'   \item{predictions}{Data frame with sample IDs, UMAP coordinates, true labels, predicted classes, and thresholded assignments}
#'   \item{probability_threshold}{Probability threshold used for "Other" assignment}
#'
#' @examples
#' result <- DLBCLone_train_mixture_model(umap_out)
#' head(result$predictions)
#' @import mclust
#'
#' @export
DLBCLone_train_mixture_model = function(umap_out,
                                        probability_threshold = 0.5,
                                        density_max_threshold = 0.05,
                                        truth_column = "lymphgen",
                                        cohort = NULL,
                                        truth_classes = c("EZB","MCD","ST2","N1","BN2","Other")
                                        ){
  
  df  = umap_out$df %>% select(sample_id,!!sym(truth_column),V1,V2) %>% filter(!!sym(truth_column) != "Other")
  df = filter(df,!!sym(truth_column) %in% truth_classes)
  #mont_test_proj = montreal_gambl_c_mu$df %>% select(sample_id,lymphgen,V1,V2)

  df[[truth_column]] <- as.factor(df[[truth_column]])


  #model without Others
  gmm_supervised <- MclustDA(
    df[, c("V1", "V2")],
    class = df[[truth_column]],
    modelType = "MclustDA"   # instead of "EDDA"
  )
df  = umap_out$df %>% select(sample_id,!!sym(truth_column),V1,V2) #%>% filter(!!sym(truth_column) != "Other")

pred <- predict(gmm_supervised, newdata = df[, c("V1", "V2")])


densities_list <- lapply(names(gmm_supervised$models), function(class_name) {
  model <- gmm_supervised$models[[class_name]]
  
  modelName <- model$modelName
  params <- model$parameters
  
  # Compute log-densities for each component
  logdens <- cdens(
    data = df[, c("V1", "V2")],
    modelName = modelName,
    parameters = params,
    logarithm = TRUE
  )
  
  # cdens can return a matrix (samples x components); sum components
  class_density <- exp(logdens)           # back to raw densities
  if (is.matrix(class_density)) {
    class_density <- rowSums(class_density)
  }
  
  class_density
})

densities <- do.call(cbind, densities_list)
colnames(densities) <- names(gmm_supervised$models)

max_density <- apply(densities, 1, max)
max_prob <- apply(pred$z, 1, max)


df$DLBCLone_g = pred$classification


df$DLBCLone_go <- ifelse(
  (max_prob < probability_threshold) | (max_density < density_max_threshold),
  "Other",
  as.character(pred$classification)
)


df$cohort = cohort
to_return = list(gaussian_mixture_model = gmm_supervised,
                 predictions = df,
                 truth_classes = truth_classes,
                 truth_column = truth_column,
                 probability_threshold = probability_threshold
                 )

return(to_return)
}

#' Predict DLBCLone Class Membership Using a Trained Gaussian Mixture Model
#'
#' Applies a previously trained supervised Gaussian mixture model (GMM) to UMAP-projected data for DLBCLone subtypes.
#' Assigns class predictions and optionally reclassifies samples as "Other" based on probability and density thresholds.
#'
#' @param model Fitted \code{MclustDA} model object, as returned by \code{DLBCLone_train_mixture_model}.
#' @param umap_out List. Output from \code{make_and_annotate_umap}, containing a data frame with UMAP
#' coordinates for the samples to be classified with the model. This must be projected using the same UMAP model
#' that was generated using the training data.
#' @param probability_threshold Numeric. Minimum posterior probability required to assign a class (default: 0.5).
#' @param density_max_threshold Numeric. Minimum maximum density required to assign a class (default: 0.05).
#' @param cohort Optional character. Cohort label to annotate predictions.
#'
#' @details
#' - Uses the provided \code{MclustDA} model to predict class membership for each sample in the UMAP projection.
#' - Computes per-class densities and posterior probabilities for each sample.
#' - Samples with low maximum probability or density are reclassified as "Other".
#' - Returns both raw and thresholded class assignments, respectively under the columns DLBCLone_g and DLBCLone_go.
#'
#' @return A list with:
#'   \item{gaussian_mixture_model}{Fitted \code{MclustDA} model object}
#'   \item{predictions}{Data frame with sample IDs, UMAP coordinates, predicted classes, and thresholded assignments}
#'   \item{probability_threshold}{Probability threshold used for "Other" assignment}
#'
#' @examples
#' # Predict on new UMAP data using a trained mixture model:
#' result <- DLBCLone_predict_mixture_model(model, umap_out)
#' head(result$predictions)
#'
#' @import mclust
#' @export
DLBCLone_predict_mixture_model = function(model,
                                          umap_out,
                                        probability_threshold = 0.5,
                                        density_max_threshold = 0.05,
                                        cohort = NULL
                                        ){
  df  = umap_out$df %>% select(sample_id,V1,V2) 
  gmm_supervised = model
pred <- predict(gmm_supervised, newdata = df[, c("V1", "V2")])

densities_list <- lapply(names(gmm_supervised$models), function(class_name) {
  model <- gmm_supervised$models[[class_name]]
  modelName <- model$modelName
  params <- model$parameters
  
  # Compute log-densities for each component
  logdens <- cdens(
    data = df[, c("V1", "V2")],
    modelName = modelName,
    parameters = params,
    logarithm = TRUE
  )
  
  # cdens can return a matrix (samples x components); sum components
  class_density <- exp(logdens)           # back to raw densities
  if (is.matrix(class_density)) {
    class_density <- rowSums(class_density)
  }
  
  class_density
})

densities <- do.call(cbind, densities_list)
colnames(densities) <- names(gmm_supervised$models)

max_density <- apply(densities, 1, max)
max_prob <- apply(pred$z, 1, max)


df$DLBCLone_g = pred$classification


df$DLBCLone_go <- ifelse(
  (max_prob < probability_threshold) | (max_density < density_max_threshold),
  "Other",
  as.character(pred$classification)
)


df$cohort = cohort
to_return = list(gaussian_mixture_model = gmm_supervised,
                 predictions = df,
                 probability_threshold = probability_threshold
                 )
return(to_return)
}


# Some comment here
# 

old_process_votes <- function(df,
                          raw_col = "label",
                          group_labels = c("EZB", "MCD", "ST2", "BN2", "N1", "Other"),
                          vote_labels_col = "vote_labels",
                          k,
                          other_vote_multiplier = 2,
                          score_purity_requirement = 1,
                          weighted_votes_col = "weighted_votes") {
  if(missing(k)){
    stop("k value is required")
  }
  score_thresh = 2 * k

  count_labels_in_string <- function(string, labels) {
    tokens <- str_split(string, ",")[[1]]
    map_int(labels, ~ sum(tokens == .x))
  }

  extract_weighted_scores <- function(label_str, vote_str, labels) {
    lbls  <- str_split(label_str, ",")[[1]]
    votes <- as.numeric(str_split(vote_str, ",")[[1]])
    map_dbl(labels, ~ sum(votes[lbls == .x])) %>%
      set_names(paste0(labels, "_score"))
  }

  get_top_score_group <- function(label_str, vote_str, labels) {
    lbls  <- str_split(label_str, ",")[[1]]
    votes <- as.numeric(str_split(vote_str, ",")[[1]])
    scores_by_label <- set_names(map_dbl(labels, ~ sum(votes[lbls == .x])), labels)
    top    <- names(scores_by_label)[which.max(scores_by_label)]
    value  <- scores_by_label[[top]]
    list(top_score_group = top, top_group_score = value)
  }

  df_out <- df %>%
    mutate(.id = row_number()) %>%
    rowwise() %>%
    mutate(
      # 1) counts
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

      # 2) scores
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

      # 3) top score group summary
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
    ungroup() %>%
    unnest_wider(counts) %>%
    unnest_wider(scores) %>%
    unnest_wider(score_summary) %>%
    rowwise() %>%
mutate(
  top_group_count = get(paste0(top_group, "_NN_count"))
) %>%
ungroup()

  # Optional adjustments based on external columns (if present)
  if ("neighbors_other" %in% colnames(df)) {
    df_out <- df_out %>%
      mutate(Other_count = neighbors_other)
  }

  if ("other_weighted_votes" %in% colnames(df)) {
    df_out <- df_out %>%
      mutate(Other_score = other_weighted_votes)
  }

  # Optional override of top group if Other dominates by a fold
  if (all(c("top_group_count", "Other_count") %in% colnames(df_out))) {
    df_out <- df_out %>%
      mutate(by_vote = top_group_count) %>%
      mutate(by_vote_opt = ifelse(top_group_count * other_vote_multiplier > Other_count, top_group, "Other"))
  }
  
  df_out = mutate(df_out,
                  by_score = top_score_group,
                  score_ratio = top_group_score / Other_score,
                  by_score_opt=ifelse(score_ratio > score_purity_requirement | top_group_score > score_thresh,top_score_group,"Other"))
  
  return(df_out)
}

