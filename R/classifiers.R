#' Classify FL samples into cFL/dFL subgroups.
#'
#' Use the random forest prediction model to assemble the binary matrix and use
#' it to classify FL tummors into cFL/dFL
#'
#' @param these_samples_metadata The metadata data frame that contains sample_id
#'      column with ids for the samples to be classified. Required argument.
#' @param maf_data The MAF data frame to be used for matrix assembling. At least
#'      must contain the first 45 columns of standard MAF format. Required
#'      argument. The maf data must be in the grch37 projection.
#' @param matrix Optionally, if the binary feature matrix is already prepared,
#'      it can be provided in this argument.
#' @param output The output to be returned after prediction is done. Can be one
#'      of predictions, matrix, or both. Defaults to predictions.
#' @param model The RF model. Classifier from the paper describing cFL is used.
#'      It is not recommended to change the value of this parameter. This is a
#'      developer-only option and used for testing and functionality improvement.
#'
#' @return data frame, binary matrix, or both
#' @export
#' @rawNamespace import(randomForest, except = c("combine"))
#' @import dplyr readr GAMBLR.data tidyr tibble GAMBLR.helpers
#'
#' @examples
#' meta <- GAMBLR.data::sample_data$meta %>%
#'     filter(pathology == "FL")
#'
#' maf <- GAMBLR.data::get_coding_ssm(
#'     these_samples_metadata = meta
#' )
#'
#' classify_fl(
#'     these_samples_metadata = meta,
#'     maf_data = maf
#' )
#' classify_fl(
#'     these_samples_metadata = meta,
#'     maf_data = maf,
#'     output = "both"
#' )
#'
classify_fl <- function(
    these_samples_metadata,
    maf_data,
    matrix,
    output = "predictions",
    model = RFmodel_FL
) {

    if(missing(matrix)){
        if(missing(these_samples_metadata) & missing(maf_data)){
            stop(
                "Exiting. Provide the sample metadata and maf data."
            )
        }

        # Ensure maf is in grch37 genome build
        maf_build = unique(maf_data$NCBI_Build)
        if(maf_build != "GRCh37"){
            stop("The provided maf data must be in the grch37 projection.")
        }
        # Ensure maf data has correct formatting
        maf_data <- maf_data %>%
            dplyr::mutate(
                Start_Position = as.numeric(Start_Position),
                End_Position = as.numeric(End_Position)
            )

        # Establish minimum required set of genes
        req_features <- rownames(
            model$importance
        )

        # SSM
        pattern <- c("HOTSPOT|_TSS|inKATdomain|_intronic|_intron_1")
        ssm_features <- req_features[!grepl(pattern, req_features)]
        # Handle HLA gene names
        ssm_features = gsub(
            "_",
            "-",
            ssm_features
        )

        # Hotspots
        hotspot_features <- req_features[grepl("HOTSPOT|inKATdomain", req_features)]

        # aSHM
        pattern <- c("_TSS|_intronic|_intron_1")
        ashm_features <- req_features[grepl(pattern, req_features)]
        ashm_features_bed <- GAMBLR.data::somatic_hypermutation_locations_GRCh37_v0.2 %>%
            dplyr::mutate(name = paste(gene, region, sep = "_")) %>%
            dplyr::mutate(
                name = case_when(
                    name == "PAX5_intron-1" ~ "PAX5_intron_1",
                    name == "PAX5_TSS-1" ~ "PAX5_TSS",
                    name == "SGK1_TSS-1" ~ "SGK1_TSS",
                    TRUE ~ name
                )
            ) %>%
            dplyr::filter(name %in% ashm_features) %>%
            dplyr::mutate(
                chrom = chr_name,
                start = hg19_start,
                end = hg19_end
            ) %>%
            dplyr::select(chrom, start, end, name) %>%
            dplyr::mutate(
                chrom = gsub("chr", "", chrom)
            )


        # Simplify gene names for all features to tabulate their mut status
        req_features <- gsub(
            "HOTSPOT|_TSS|inKATdomain|_intronic|_intron_1", # drop any extra info
            "",
            req_features
        )
        req_features <- gsub(
            "_", # HLA gene name handling
            "-",
            req_features
        )
        req_features <- sort(
            unique(
                req_features
            )
        )

        # Handle possible missing samples in maf compared to metadata
        found_samples <- length(unique(maf_data$Tumor_Sample_Barcode))
        requested_samples <- length(unique(these_samples_metadata$sample_id))
        if(!found_samples == requested_samples){
            missing_samples <- setdiff(
                    unique(these_samples_metadata$sample_id),
                    unique(maf_data$Tumor_Sample_Barcode)
                )
            message(
                paste0(
                    "Did not find SSM for all samples. Only the data for ",
                    found_samples,
                    " was available in the maf. The missing samples are: "
                )
            )
            message(
                paste(
                    as.character(missing_samples),
                    collapse = ", "
                )
            )
            message(
                "They will be dropped since their mutation data is not provided."
            )
            # Drop missing samples from metadata
            these_samples_metadata <- these_samples_metadata %>%
                dplyr::filter(
                    sample_id %in% maf_data$Tumor_Sample_Barcode
                )
        }else{
            message(
                "The SSM for all samples were found in the provided data."
            )
        }
        message(
            "Proceeding to matrix assembling..."
        )

        # Generate binary matrix for SSMs and hotspots
        ssm_matrix <- tabulate_ssm_status(
            gene_symbols = ssm_features,
            these_samples_metadata = these_samples_metadata,
            maf_data = maf_data,
            genome_build = "grch37",
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
                    "CREBBPinKATdomain" = "CREBBPHOTSPOT"
                )
        }


        # Generate binary matrix for ashm
        overlap <- cool_overlaps(
            data1 = maf_data,
            data2 = ashm_features_bed,
            columns2 = c("chrom", "start", "end")
        )

        ashm_matrix <- overlap %>%
            dplyr::group_by(Tumor_Sample_Barcode, name) %>%
            dplyr::tally() %>%
            dplyr::ungroup() %>%
            tidyr::pivot_wider(
                names_from = name,
                values_from = n
            ) %>%
            tibble::column_to_rownames("Tumor_Sample_Barcode") %>%
            replace(is.na(.), 0)

        # Handle possible samples with no muts
        ashm_matrix <- complete_missing_from_matrix(
            ashm_matrix,
            these_samples_metadata$sample_id
        )


        ashm_matrix[ashm_matrix <= 5] = 0
        ashm_matrix[ashm_matrix > 5] = 1

        ashm_matrix <- check_for_missing_features(
            ashm_matrix,
            ashm_features
        )


        # Combine together the SSM, hotspots, and aSHM
        assembled_matrix <- rownames_to_column(ssm_matrix, "sample_id") %>%
            left_join(ashm_matrix %>% rownames_to_column("sample_id"), by = "sample_id") %>%
            column_to_rownames("sample_id")

            # Check for missing features
        assembled_matrix <- check_for_missing_features(
            assembled_matrix,
            c(
                ssm_features,
                hotspot_features,
                ashm_features
            )
        )
    }else{
        assembled_matrix <- matrix
    }

    if("HLA-DMB" %in% colnames(assembled_matrix)){
        assembled_matrix <- assembled_matrix %>%
            dplyr::rename(
                "HLA_DMB" = "HLA-DMB"
            )
    }

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
        type = "Vote"
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

    if(output == "predictions"){
        return(prediction)
    }else if (output == "matrix") {
       return(assembled_matrix)
    }else if(output == "both"){
        return(
            list(
                predictions = prediction,
                matrix = assembled_matrix)
        )
    }else{
        stop(
            "Invalid output type. Please specify predictions, matrix, or both."
        )
    }

}

#' Classify DLBCLs according to genetic subgroups.
#'
#' Using the user-provided or GAMBLR.data-retrieved data, this function will assemble the matrix according to the approach of
#' Chapuy et al (2018) or Lacy et al (2020) classifiers. Since neither of this classifiers is publicly released, we have implemented a solution
#' that closely (> 92% accuracy) recapitulates each of these systems. For the classifier of Chapuy et al, the constructed matrix will be
#' used to calculate class probability using the bundled feature weights obtained from our reproduction of the classifier. For the Lacy et al
#' classifier, the matrix will be used for prediction of random forest model, which is supplied with the GAMBLR.predict package. Following the modification
#' of Lacy classifier described in Runge et al (PMID 33010029), specifying the method of this function as hmrn will also consider
#' truncating mutations in NOTCH1 for the separate N1 subgroup.
#'
#' @param these_samples_metadata The metadata data frame that contains sample_id column with ids for the samples to be classified.
#' @param maf_data The MAF data frame to be used for matrix assembling. At least must contain the first 45 columns of standard MAF format.
#' @param only_maf_data Whether to restrict matrix generation to maf data only. Only supported in the lymphgenerator mode. Default is FALSE (use SV and CNV).
#' @param seg_data The SEG data frame to be used for matrix assembling. Must be of standard SEG formatting, for example, as returned by get_sample_cn_segments.
#' @param sv_data The SV data frame to be used for matrix assembling. Must be of standard BEDPE formatting, for example, as returned by get_combined_sv.
#' @param this_seq_type The seq_type of the samples. Only used to retrieve data through GAMBLR.data when it is not provided. Defaults to genome.
#' @param projection The projection of the samples. Only used to retrerive data through GAMBLR.data when it is not provided. Defaults to grch37.
#' @param output The output to be returned after the prediction is done. Can be one of predictions, matrix, or both. Defaults to both.
#' @param method Classification method. One of chapuy (used as default), lacy, or hmrn.
#' @param adjust_ploidy Whether to perform ploidy adjustment for the CNV data. Defaults to TRUE (recommended).
#' @param annotate_sv Whether to perform SV annotation on the supplied SV data frame. Defaults to TRUE.
#'
#' @return data frame with classification, binary matrix used in classification, or both
#' @export
#' @import dplyr readr
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
    only_maf_data = FALSE,
    seg_data,
    sv_data,
    this_seq_type = "genome",
    projection = "grch37",
    output = "both",
    method = "chapuy",
    adjust_ploidy = TRUE,
    annotate_sv = TRUE
){
    # If no metadata is provided, just get all DLBCLs
    if(missing(these_samples_metadata)){
        message("No metadata is provided.")
        stop(
            "Please provide sample metadata."
        )
    }else{
        these_samples_metadata <- these_samples_metadata %>%
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
    if(!only_maf_data & missing(seg_data)){
        message("No CNV data is provided.")
        message("Will retreive segments available through GAMBL.")

        seg_data = get_sample_cn_segments(
            these_samples_metadata = these_samples_metadata,
            projection = projection)
    }

    # If no SV data is provided, get the SVs from GAMBL
    if(!only_maf_data & missing(sv_data) & method %in% c("chapuy", "lymphgenerator")){
        message("No SV data is provided.")
        message("Will retreive SVs available through GAMBL.")

        sv_data <- get_manta_sv() %>%
                dplyr::filter(
                tumour_sample_id %in% these_samples_metadata$sample_id
            )
    }

    if(!only_maf_data){
        if(adjust_ploidy){
            seg_data <- adjust_ploidy(
                seg_data %>% rename("sample"="ID"),
                projection = projection
            )
        }

        seg_data <- seg_data %>%
            as.data.frame

        if(annotate_sv){
            sv_data <- sv_data %>%
                annotate_sv(
                    genome_build = projection
                ) %>%
                dplyr::filter(!is.na(partner))
        }
    }

    # Now use this data to classify the samples according to one of the systems
    if(method == "chapuy"){

      predictions <- classify_dlbcl_chapuy(
          these_samples_metadata,
          maf_data,
          seg_data,
          sv_data,
          projection = projection,
          output = output
      )

    }else if (method == "lacy") {
      predictions <- classify_dlbcl_lacy(
          these_samples_metadata,
          maf_data,
          seg_data,
          projection = projection,
          output = output
      )

    }else if (method == "hmrn") {
      predictions <- classify_dlbcl_lacy(
          these_samples_metadata,
          maf_data,
          seg_data,
          projection = projection,
          output = output,
          include_N1 = TRUE
      )
    }else if (method == "lymphgenerator") {
        if(only_maf_data){
            predictions <- classify_dlbcl_lymphgenerator(
                these_samples_metadata,
                maf_data,
                projection = projection,
                output = output,
                seq_type = this_seq_type
            )
        }else{
            predictions <- classify_dlbcl_lymphgenerator(
                these_samples_metadata,
                maf_data,
                sv_data,
                seg_data,
                projection = projection,
                output = output,
                seq_type = this_seq_type
            )
        }
    }

    else{
      stop("You requested an unvalid method. Please choose one of chapuy, lacy or hmrn")
    }

    return(predictions)

}

#' Classify BL samples into genetic subgroups.
#'
#' Use the random forest prediction model to assemble the binary matrix and use it to classify BL tummors into genetic subgroups
#'
#' @param these_samples_metadata The metadata data frame that contains sample_id column with ids for the samples to be classified.
#' @param maf_data The MAF data frame to be used for matrix assembling. At least must contain the first 45 columns of standard MAF format.
#' @param this_seq_type The seq_type of the samples. Only really used to retrerive mutations when maf data is not provided and to be retreived through GAMBLR.data. Defaults to genome.
#' @param projection The projection of the samples. Defaults to grch37.
#' @param output The output to be returned after prediction is done. Can be one of predictions, matrix, or both. Defaults to both.
#' @param ashm_cutoff Numeric value indicating number of mutations for binarizing aSHM feature. Recommended to use the default value (3).
#'
#' @return data frame with classification, binary matrix used in classification, or both
#' @export
#' @rawNamespace import(randomForest, except = c("combine"))
#' @import dplyr GAMBLR.data tidyr tibble
#'
#' @examples
#' test_meta <- get_gambl_metadata(case_set = "BL-DLBCL-manuscript-HTMCP")
#' predictions <- classify_bl(these_samples_metadata = test_meta)
#' predictions <- classify_bl(these_samples_metadata = test_meta, output = "predictions")
#'
classify_bl <- function(
    these_samples_metadata,
    maf_data,
    this_seq_type = "genome",
    projection = "grch37",
    output = "both",
    ashm_cutoff = 3
){
    # If no metadata is provided, just get all BLs
    if(missing(these_samples_metadata)){
        message("No metadata is provided.")
        stop(
            "Please provide sample metadata."
        )
    }else{
        these_samples_metadata <- these_samples_metadata %>%
            dplyr::filter(pathology == "BL")
    }


    projection <- handle_genome_build(projection)

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

    # Now use this data to classify the samples according to one of the systems
    # Generate binary matrix for SSMs
    ssm_matrix <- get_coding_ssm_status(
        gene_symbols = c(
            rownames(RFmodel_BL$importance),
            "HLA-A", "HLA-B", "HLA-DMB"
        ),
        these_samples_metadata = these_samples_metadata,
        maf_data = maf_data,
        include_hotspots = FALSE
    )

    ssm_matrix <- ssm_matrix %>%
        column_to_rownames("sample_id")

    # Generate binary matrix for aSHM
    ashm_bed <- base::get(
        paste0(
            projection, "_ashm_regions"
        )
    ) %>%
    dplyr::filter(name %in% c(
        "LPP_TSS-1",
        "LPP-TSS-1"
        )
    )

    # Generate binary matrix for ashm
    ashm_matrix <- get_ashm_count_matrix(
        maf_data = maf,
        these_samples_metadata = these_samples_metadata,
        regions_bed = ashm_bed
    )

    colnames(ashm_matrix) <- c("LPPTSS1")

    ashm_matrix <- ashm_matrix %>%
        mutate("LPPTSS1" = ifelse(
            LPPTSS1 > ashm_cutoff,
            1,
            0
        )
    )

    # Combine together the SSM, hotspots, and aSHM
    assembled_matrix <- bind_cols(
            ssm_matrix,
            ashm_matrix
        )

    # Prepare for random forest
    colnames(assembled_matrix) <- (gsub(
            '_|-|\\.',
            '',
            colnames(assembled_matrix)
        )
    )

    # Check for missing features
    assembled_matrix <- check_for_missing_features(
        assembled_matrix,
        rownames(RFmodel_BL$importance)
    )

    # Ensure consistent ordering
    assembled_matrix <- assembled_matrix %>%
        select(
            rownames(
                RFmodel_BL$importance
            )
        )

    # Make prediction
    prediction <- predict(
        RFmodel_BL,
        assembled_matrix,
        type = "Vote"
    )

    prediction <- bind_cols(
        prediction,
        predict(
                RFmodel_BL,
                assembled_matrix
            ) %>%
            as.data.frame %>%
            `names<-`("BL_subgroup")
        ) %>%
    rownames_to_column("sample_id")

    if(output == "predictions"){
        return(prediction)
    }else if (output == "matrix") {
       return(assembled_matrix)
    }else if(output == "both"){
        return(
            list(
                predictions = prediction,
                matrix = assembled_matrix)
        )
    }else{
        stop("Invalid output type. Please specify predictions, matrix, or both.")
    }

}


#' Complete samples missing from matrix.
#'
#' If some samples are missing from the matrix, add them with filled in 0 as value and normalize their ordering for consistency.
#'
#' @param incoming_matrix A matrix or data frame that should be filled. Required parameter.
#' @param list_of_samples Vector specifying all desired samples to be present in the resulting matrix. Required parameter.
#' @param fill_in_values Value that will be used to fill in the matrix.
#' @param normalize_order Logical parameter specifying whether sample order should be according to the supplied list. Default is TRUE.
#' @param samples_in_rows Logical argument indicating whether samples are in rows or columns. Default assumes samples are in rows and columns are features.
#'
#' @return A data frame with maintained orientation (rows and columns) where samples from the supplied list are present and reordered according to the specified order.
#' @export
#' @import dplyr
#'
#' @examples
#' partial_matrix = get_coding_ssm_status(these_samples_metadata = (get_gambl_metadata(case_set = "BL--DLBCL") %>% filter(pairing_status == "unmatched")), include_hotspots = FALSE)
#' complete_matrix = complete_missing_from_matrix(partial_matrix, get_gambl_metadata() %>% pull(sample_id))
#'
complete_missing_from_matrix = function(
    incoming_matrix,
    list_of_samples,
    fill_in_values = 0,
    normalize_order = TRUE,
    samples_in_rows = TRUE
){

    # check for required arguments
    if (missing(incoming_matrix)){
        stop("Please provide initial matrix to fill.")
    }

    if (missing(list_of_samples)){
        stop("Please provide list of samples to complete the matrix and normalize order.")
    }

    # is samples are in columns, transpose the matrix so code below is generalizable
    if(!samples_in_rows){
        incoming_matrix = as.data.frame(incoming_matrix) %>%
        t
    }

    matrix_with_all_samples = rbind(
        incoming_matrix,
        matrix(
            fill_in_values:fill_in_values, # populate matrix with all 0
            length(
                setdiff(
                    list_of_samples,
                    rownames(incoming_matrix)
                )
            ), # how many rows
            ncol(incoming_matrix), # how many columns
            dimnames = list(
                setdiff(
                    list_of_samples,
                    rownames(incoming_matrix)
                ), # name rows with sample IDs
                colnames(incoming_matrix) # name columns with gene names
            )
        ) %>%
        as.data.frame
    )

    # this is very helpful in clustering
    if(normalize_order){
        matrix_with_all_samples =
            matrix_with_all_samples[
                order(match(
                    rownames(matrix_with_all_samples),
                    list_of_samples
                    )
                )
                ,
            ]
    }

    # transpose matrix back to the initial format
    # supplied by user (samples in columns)
    if(!samples_in_rows){
        matrix_with_all_samples = as.data.frame(matrix_with_all_samples) %>%
        t()
    }
    return(matrix_with_all_samples)
}


#' Will prepare the data frame of binary matrix to be used as NMF input. This means that for the features with SSM and CNV,
#' they will be squished together as one feature named GeneName-MUTorAMP or GeneName-MUTorLOSS, so the CNV features in the input data frame are expected
#' to be named GeneName_AMP or GeneName_LOSS. Next, for the genes with hotspot mutations labelled in the input data as
#' GeneNameHOTSPOT, the feature for hotspot mutation will be given preference and SSM with/without CNV will be set to 0 for that sample.
#' The naming scheme of the features as in this description is important, because the function uses regex to searh for these patters as specified.
#' Finally, if any features are provided to be dropped explicitly, they will be removed, and then the features not meeting the specified minimal
#' frequency will be removed, as well as any samples with 0 features.
#' Consistent with NMF input, in the input data frame each row is a feature, and each column is a sample. The input is expected to be numeric 1/0 with row and column names.
#'
#'
#' @param incoming_data Input data frame or matrix to prepare for NMF.
#' @param blacklisted_cnv_regex Regular expression to match in feature names when considering SSM/CNV overlap.
#' @param drop_these_features Optional argument with features to drop from resulting matrix.
#' @param min_feature_percent Minimum frequency for the feature to be returned in the resulting matrix. By default, features present in less than 0.5% of samples will be discarded.
#'
#' @return A matrix compatible with NMF input.
#' @export
#' @import dplyr
#'
#' @examples
#' data = system.file("extdata", "sample_matrix.tsv", package = "GAMBLR.predict") %>% read_tsv() %>% column_to_rownames("Feature")
#' NMF_input = massage_matrix_for_clustering(data)
#'

massage_matrix_for_clustering = function(
    incoming_data,
    blacklisted_cnv_regex = "3UTR|SV|HOTSPOT|TP53BP1|intronic",
    drop_these_features,
    min_feature_percent = 0.005
){

    # if there is a CNV and mutation at the same gene, squish these features together
    message("Searching for overlapping CNV and mutation features to squish together ...")
    feat_with_cnv_data = rownames(incoming_data)[grepl("AMP|LOSS", rownames(incoming_data))]
    output_data = incoming_data

    for (g in feat_with_cnv_data){
        this_feature = unlist(strsplit(g, split='_', fixed=TRUE))
        red_features <- rownames(output_data)[grepl(this_feature[1], rownames(output_data))]
        red_features <- red_features[!grepl(blacklisted_cnv_regex, red_features)] # these features to be kept separately
        if(length(red_features)>1){
        message(paste0("Found redundant features for gene ", red_features[1], ", processing ..."))
        output_data[,output_data[red_features[2],]>0][red_features,][red_features[1],] = 1
        rownames(output_data)[rownames(output_data)==red_features[1]] = paste0(this_feature[1],
                                                                "-MUTor",
                                                                this_feature[2])
        output_data = output_data[!rownames(output_data) %in% red_features[2],]

        }
    }

    message("Success")

    # if there is a hotspot and SSM for same gene, give priority to hotspot
    message("Searching for overlapping HOTSPOT and mutation features to squish together ...")
    feat_with_hotspot_data = rownames(output_data)[grepl("HOTSPOT", rownames(output_data))]
    for (hot in feat_with_hotspot_data){
        this_gene=gsub("HOTSPOT","", hot)
        # this gene may also have CNV data already squished
        maybe_cnv = grepl("MUTor",
                        rownames(output_data[grepl(this_gene,
                                            rownames(output_data)),]))
        if("TRUE" %in% maybe_cnv){ # if it has the cnv data then use the name of gene with LOSS or AMP respectively
        this_gene = rownames(output_data[grepl(this_gene, rownames(output_data)),])[maybe_cnv]
        message(paste0("Found hotspot for gene ", this_gene, " that also has CNV data, processing ..."))
        output_data[,(output_data[c(this_gene),]>0 & output_data[c(hot),]==1)][c(this_gene, hot),][c(this_gene),] = 0
        }else{ # otherwise just use the gene nae
        message(paste0("Found hotspot for gene ", this_gene, ", processing ..."))
        output_data[,(output_data[c(this_gene),]>0 & output_data[c(hot),]==1)][c(this_gene, hot),][c(this_gene),] = 0
        }
        # if the above statement work, then there should be no overlaps between hotspot and any other mutations
        # for the same gene
        if(length(output_data[,(output_data[c(this_gene),]>0 & output_data[c(hot),]==1)][c(this_gene, hot),][c(this_gene),])==0){
        message("Success")
        }else{
        message(paste0("Problem occured with the ", feat_with_hotspot_data, " and the gene ", this_gene, "and there is still overlap between mutation and hotspot."))
        break
        }
    }

    # did user provide any features they would like to drop from matrix?
    if(!missing(drop_these_features)){
        message("You provided features to be dropped from matrix, removing them ...")
        output_data =
        output_data[!rownames(output_data) %in% drop_these_features,]
        message("Success")
    }

    # drop features that are occuring at a very low frequency
    low_feat = which(rowSums(output_data) <= floor(ncol(output_data)*min_feature_percent))
    if (length(low_feat)>0){
        message(paste0 ("There are ", length(low_feat), " features underrepresented and not meeting the minimum frequency of ", min_feature_percent))
        print(names(low_feat))
        output_data = output_data[-c(low_feat),]
    }else{
        message(paste0 ("There are ", length(low_feat), " features not meeting the minimum frequency of ", min_feature_percent))
        message("Proceeding without dropping any feature ...")
    }

    # are there any samples with 0 features? Yes, 1 exome and 1 genome
    samples_with_zero_feat = which(colSums(output_data) == 0)
    if (length(samples_with_zero_feat)>0){
        message(paste0 ("There are ", length(samples_with_zero_feat), " samples with no features and they will be dropped from matrix: "))
        print(names(samples_with_zero_feat))
        output_data = output_data[, -c(samples_with_zero_feat)]
    }else{
        message("All samples in your matrix are having at least one feature. Proceeding without dropping any samples ...")
    }

    # convert to matrix explicitly to make it NMF-input compatible
    output_data = as.matrix(output_data)

    return(output_data)

}
