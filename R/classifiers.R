#' Classify FL samples into cFL/dFL subgroups.
#'
#' Use the random forest prediction model to assemble the binary matrix and use it to classify FL tummors into cFL/dFL
#'
#' @param these_samples_metadata The metadata data frame that contains sample_id column with ids for the samples to be classified.
#' @param maf_data The MAF data frame to be used for matrix assembling. At least must contain the first 45 columns of standard MAF format.
#' @param model The RF model. Classifier from the paper describing cFL is used. It is not recommended to change the value of this parameter.
#' @param this_seq_type The seq_type of the samples. Only really used to retrerive mutations when maf data is not provided and to be retreived through GAMBLR. Defaults to genome.
#' @param output The output to be returned after prediction is done. Can be one of predictions, matrix, or both. Defaults to predictions.
#'
#' @return data frame with classification, binary matrix used in classification, or both
#' @export
#' @import dplyr readr stringr randomForest GAMBLR
#'
#' @examples
#' test_meta <- get_gambl_metadata(case_set="tFL-study")
#' predictions = classify_fl(these_samples_metadata = test_meta)
#' predictions = classify_fl(these_samples_metadata = test_meta, output = "both")
#'
classify_fl <- function(
    these_samples_metadata,
    maf_data,
    model = RFmodel_FL,
    this_seq_type = "genome",
    output = "predictions"
) {

    # Establish minimum required set of genes
    req_features <- rownames(
        model$importance
    )

    ssm_features <- req_features[!grepl("HOTSPOT|_TSS|inKATdomain|_intronic|_intron_1", req_features)]
    ssm_features = gsub(
        "_",
        "-",
        ssm_features
    )
    hotspot_features <- req_features[grepl("HOTSPOT|inKATdomain", req_features)]
    ashm_features <- req_features[grepl("_TSS|_intronic|_intron_1", req_features)]
    ashm_features_bed <- grch37_ashm_regions %>%
        dplyr::mutate(
            name = gsub("-", "_", grch37_ashm_regions$name)
        ) %>%
        dplyr::filter(name %in% ashm_features | name %in% c("PAX5_TSS_1", "SGK1_TSS_1"))

    req_features <- gsub(
        "HOTSPOT|_TSS|inKATdomain|_intronic|_intron_1",
        "",
        req_features
    )
    req_features <- gsub(
        "_",
        "-",
        req_features
    )
    req_features <- sort(
        unique(
            req_features
        )
    )

    if(missing(these_samples_metadata) & missing(maf_data)){
        stop("Exiting. Please provide the sample metadata or maf data to use in classification.")
    }else if (missing(maf_data)){
       message(
            "No maf data was provided. Retreiving SSMs using GAMBLR..."
       )
       maf_data =  get_ssm_by_samples(
            these_samples_metadata = these_samples_metadata,
            seq_type = this_seq_type,
            subset_from_merge = TRUE,
            augmented = FALSE
       )
       found_samples <- length(unique(maf_data$Tumor_Sample_Barcode))
       requested_samples <- length(unique(these_samples_metadata$sample_id))
       if(!found_samples == requested_samples){
            message(
                paste0(
                    "Did not find SSM for all samples. Only the data for ",
                    found_samples,
                    " was available through GAMBLR. The missing samples are: "
                )
            )
            message(
                setdiff(
                    unique(these_samples_metadata$sample_id),
                    unique(maf_data$Tumor_Sample_Barcode)
                )
            )
            # Drop missing samples from metadata
            these_samples_metadata <- these_samples_metadata %>%
                dplyr::filter(
                    sample_id %in% maf_data$Tumor_Sample_Barcode
                )
       }else{
            message(
                "The SSM for all samples were found in GAMBLR. Proceeding to matrix assembling."
            )
       }
    }else if (missing(these_samples_metadata)){
        message(
            "The metadata was not provided. Retreiving the metadata through GAMBLR..."
        )
        these_samples_metadata <- get_gambl_metadata(seq_type_filter = this_seq_type) %>%
            filter(sample_id %in% maf_data$Tumor_Sample_Barcode)

    }else{
        message(
            "Using the provided metadata and maf to assemble the matrix..."
        )
    }

    # Generate binary matrix for SSMs and hotspots
    ssm_matrix <- get_coding_ssm_status(
        gene_symbols = ssm_features,
        these_samples_metadata = these_samples_metadata,
        maf_data = maf_data,
        genes_of_interest = gsub(
            "HOTSPOT|inKATdomain",
            "",
            hotspot_features
        )
    )

    ssm_matrix <- ssm_matrix %>%
        column_to_rownames("sample_id")

    if("CREBBPHOTSPOT" %in% colnames(ssm_matrix)){
        ssm_matrix <- ssm_matrix %>%
            dplyr::rename(
                "CREBBPinKATdomain" = "CREBBPHOTSPOT",
                "HLA_DMB" = "HLA-DMB"
            )
    }

    # Generate binary matrix for ashm
    ashm_matrix <- get_ashm_count_matrix(
        ashm_features_bed,
        maf_data = maf_data,
        these_samples_metadata = these_samples_metadata,
        seq_type = this_seq_type
    )

    ashm_matrix[ashm_matrix<=5] = 0
    ashm_matrix[ashm_matrix>5] = 1


    ashm_matrix <- ashm_matrix %>%
        rename(
            "PAX5_TSS" = "PAX5_TSS_1",
            "SGK1_TSS" = "SGK1_TSS_1"
    )


    # Combine together the SSM, hotspots, and aSHM
    assembled_matrix <- bind_cols(
            ssm_matrix,
            ashm_matrix
        )

    # Check for missing features
    assembled_matrix <- check_for_missing_features(
        assembled_matrix,
        c(
          ssm_features,
          hotspot_features,
          ashm_features
        )
    )

    # Ensure consistent ordering
    assembled_matrix <- assembled_matrix %>%
        select(
            rownames(
                model$importance
            )
        )

    # Make prediction
    prediction <- predict(
        model,
        assembled_matrix,
        type="Vote"
    )

    prediction <- bind_cols(
        prediction,
        predict(
                model,
                assembled_matrix
            ) %>%
            as.data.frame %>%
            `names<-`("is_cFL")
        ) %>%
    rownames_to_column("sample_id")

    if(output=="predictions"){
        return(prediction)
    }else if (output=="matrix") {
       return(assembled_matrix)
    }else if(output=="both"){
        return(
            list(
                predictions = prediction,
                matrix = assembled_matrix)
        )
    }else{
        stop("Invalid output type. Please specify predictions, matrix, or both.")
    }


}

#' Classify DLBCLs according to genetic subgroups.
#'
#' Using the user-provided or GAMBLR-retrieved data, this function will assemble the matrix according to the approach of
#' Chapuy et al (2018) or Lacy et al (2020) classifiers. Since neither of this classifiers is publicly released, we have implemented a solution
#' that closely (> 92% accuracy) recapitulates each of these systems. For the classifier of Chapuy et al, the constructed matrix will be
#' used to calculate class probability using the bundled feature weights obtained from our reproduction of the classifier. For the Lacy et al
#' classifier, the matrix will be used for prediction of random forest model, which is supplied with the GAMBLR package. Following the modification
#' of Lacy classifier described in Runge et al (PMID 33010029), specifying the method of this function as hmrn will also consider
#' truncating mutations in NOTCH1 for the separate N1 subgroup.
#'
#' @param these_samples_metadata The metadata data frame that contains sample_id column with ids for the samples to be classified.
#' @param maf_data The MAF data frame to be used for matrix assembling. At least must contain the first 45 columns of standard MAF format.
#' @param seg_data The SEG data frame to be used for matrix assembling. Must be of standard SEG formatting, for example, as returned by get_sample_cn_segments.
#' @param sv_data The SV data frame to be used for matrix assembling. Must be of standard BEDPE formatting, for example, as returned by get_combined_sv.
#' @param this_seq_type The seq_type of the samples. Only used to retrerive data through GAMBLR when it is not provided. Defaults to genome.
#' @param projection The projection of the samples. Only used to retrerive data through GAMBLR when it is not provided. Defaults to grch37.
#' @param output The output to be returned after prediction is done. Can be one of predictions, matrix, or both. Defaults to both.
#' @param method Classification method. One of chapuy (used as default), lacy, or hmrn.
#' @param adjust_ploidy Whether to perform ploidy adjustment for the CNV data. Defaults to TRUE (recommended).
#'
#' @return data frame with classification, binary matrix used in classification, or both
#' @export
#' @import data.table circlize dplyr reads stringr
#'
#' @examples
#' test_meta <- get_gambl_metadata(case_set = "DLBCL-unembargoed")
#' predictions_chapuy <- classify_dlbcl(these_samples_metadata = test_meta, output = "predictions")
#' predictions_lacy <- classify_dlbcl(these_samples_metadata = test_meta, method = "lacy")
#' predictions_hmrn <- classify_dlbcl(these_samples_metadata = test_meta, method = "hmrn", output = "predictions")
#' matrix_and_predictions <- classify_dlbcl(these_samples_metadata = test_meta)
#'
classify_dlbcl <- function(
    these_samples_metadata,
    maf_data,
    seg_data,
    sv_data,
    this_seq_type = "genome",
    projection = "grch37",
    output = "both",
    method = "chapuy",
    adjust_ploidy = TRUE
){
    # If no metadata is provided, just get all DLBCLs
    if(missing(these_samples_metadata)){
        message("No metadata is provided.")
        message("Will retreive metadata for all DLBCL genomes in GAMBL.")
        these_samples_metadata <- get_gambl_metadata(
            seq_type_filter = this_seq_type
        ) %>%
        dplyr::filter(pathology == "DLBCL")
    }

    # If no maf data is provided, get the SSMs from GAMBL
    if(missing(maf_data)){
        message("No maf data is provided.")
        message("Retreiving the mutations data from GAMBL...")
        maf_data =  get_ssm_by_samples(
            these_samples_metadata = these_samples_metadata,
            seq_type = this_seq_type,
            projection = projection,
            subset_from_merge = TRUE,
            augmented = FALSE
        )
    }


    # Confirm all samples have mutations
    found_samples <- length(unique(maf_data$Tumor_Sample_Barcode))
    requested_samples <- length(unique(these_samples_metadata$sample_id))

    if(!found_samples == requested_samples){
        message(
            paste0(
                "WARNING! Did not find SSM for all samples. Only the data for ",
                found_samples,
                " was available in the maf. The missing samples are: "
            )
        )
        message(
            paste(
              setdiff(
                unique(these_samples_metadata$sample_id),
                unique(maf_data$Tumor_Sample_Barcode)
              ),
              collapse=", "
            )
        )
        # Drop missing samples from metadata
        these_samples_metadata <- these_samples_metadata %>%
            dplyr::filter(
                sample_id %in% maf_data$Tumor_Sample_Barcode
            )
    }else{
        message(
            "Success! The SSM for all samples were found in maf."
        )
    }


    # If no seg data is provided, get the CNVs from GAMBL
    if(missing(seg_data)){
        message("No CNV data is provided.")
        message("Will retreive segments available through GAMBL.")

        seg_data = get_sample_cn_segments(
            sample_list = these_samples_metadata$sample_id,
            multiple_samples = TRUE,
            projection = projection)
    }

    if(adjust_ploidy){
        seg_data <- adjust_ploidy(
            seg_data %>% rename("sample"="ID"),
            projection = projection
        )
    }

    seg_data <- seg_data %>%
        as.data.table %>%
        setkey(chrom, start, end)


    # Now use this data to classify the samples according to one of the systems
    if(method=="chapuy"){

      # SVs are needed only for Chapuy predictions
      # If no SV data is provided, get the SVs from GAMBL
      if(missing(sv_data)){
          message("No SV data is provided.")
          message("Will retreive SVs available through GAMBL.")

          sv_data <- get_manta_sv() %>%
              dplyr::filter(
                tumour_sample_id %in% these_samples_metadata$sample_id
              )
      }

      sv_data <- sv_data %>%
          annotate_sv(
              genome_build = projection
          ) %>%
          dplyr::filter(!is.na(partner))

      predictions <- classify_dlbcl_chapuy(
          these_samples_metadata,
          maf_data,
          seg_data,
          sv_data,
          projection = projection,
          output = output
      )

    }else if (method=="lacy") {
      predictions <- classify_dlbcl_lacy(
          these_samples_metadata,
          maf_data,
          seg_data,
          projection = projection,
          output = output
      )

    }else if (method=="hmrn") {
      predictions <- classify_dlbcl_lacy(
          these_samples_metadata,
          maf_data,
          seg_data,
          projection = projection,
          output = output,
          include_N1 = TRUE
      )
    }

    else{
      stop("You requested an unvalid method. Please choose one of chapuy, lacy or hmrn")
    }

    return(predictions)

}
