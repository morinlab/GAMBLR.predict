
#' Check matrix against missing features.
#'
#' Operate on the matrix supplied by user and check it for any features missing compared to the provided set.
#'
#' @param incoming_matrix The incoming matrix to be checked for any missing features.
#' @param feature_set The feature set.
#' @return matrix where column names correspond to all features from the provided set
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
#' Use the feature weights from NMF model to assemble the binary matrix and classify DLBCL tumors based on C0-C5 system of Chapuy et al
#'
#' @param these_samples_metadata The metadata data frame that contains sample_id column with ids for the samples to be classified.
#' @param maf_data The MAF data frame to be used for matrix assembling. At least must contain the first 45 columns of standard MAF format.
#' @param seg_data The SEG data frame to be used for matrix assembling. Must be of standard SEG formatting, for example, as returned by get_sample_cn_segments.
#' @param sv_data The SV data frame to be used for matrix assembling. Must be of standard BEDPE formatting, for example, as returned by get_combined_sv.
#' @param projection The projection of the samples. Only used to retrerive data through GAMBLR when it is not provided. Defaults to grch37.
#' @param output The output to be returned after prediction is done. Can be one of predictions, matrix, or both. Defaults to both.
#' @return data frame with classification, binary matrix used in classification, or both
#' @import data.table circlize dplyr readr stringr
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
            GAMBLR:::coding_class
        )) %>%
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
        arm_coordinates <- GAMBLR::chromosome_arms_grch37
        cytoband_coordinates <- circlize::read.cytoband(species = "hg19")$df %>%
            `names<-`(c("chr", "start", "end", "cytoband", "extra")) %>%
            dplyr::mutate(chr = gsub("chr", "", chr)) %>%
            dplyr::mutate(cytoband=paste0(chr,cytoband)) %>%
            dplyr::select(-extra)
    }else{
        arm_coordinates <- GAMBLR::chromosome_arms_hg38
        cytoband_coordinates <- circlize::read.cytoband(species = "hg38")$df %>%
            `names<-`(c("chr", "start", "end", "cytoband", "extra")) %>%
            dplyr::mutate(cytoband=paste0(chr,cytoband)) %>%
            dplyr::select(-extra)
    }

    # First the arm features
    cnv_features_arm <- arm_coordinates %>%
        mutate(arm = paste0(
                chromosome,
                arm)
            ) %>%
        left_join(
                chapuy_features$cnv_features_arm,
                .,
                by="arm"
            ) %>%
        as.data.table %>%
        setkey(chromosome, start, end)

    # Next, the cytoband features
    cnv_features_cytoband <- cytoband_coordinates %>%
        left_join(
                chapuy_features$cnv_features_cytoband,
                .,
                by="cytoband"
            ) %>%
        as.data.table %>%
        setkey(chr, start, end)

    cnv_arms <- foverlaps(
          seg_data,
          cnv_features_arm,
          nomatch = 0
        ) %>%
        dplyr::select(sample, arm, CNV, log.ratio) %>%
        dplyr::rename("feature"="arm")

    cnv_cytobands <-  foverlaps(
          seg_data,
          cnv_features_cytoband,
          nomatch = 0
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
        )) %>%
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

    # Layer in which cluster the sample belongs to
    # by taking the highest sum of weights
    predictions$predict <- colnames(predictions)[apply(predictions,1,which.max)]

    predictions <- predictions %>%
      rownames_to_column("sample_id")

    # Account for C0 samples, which will have all weights calculated as 0
    predictions <- predictions %>%
      rowwise() %>%
      dplyr::mutate(
        predict = ifelse(
          sum(C1:C5)==0,
          "C0",
          predict
      )) %>%
      ungroup %>%
      as.data.frame %>%
      dplyr::rename(
        "Chapuy_cluster"="predict"
      )

    if(output == "predictions"){
      return(predictions)
    }else if (output == "both") {
      return(
        list(
          chapuy_matrix = chapuy_feature_matrix$complete_matrix,
          chapuy_predictons = predictions
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
#' Use the random forest model to classify DLBCL tumors based on system of Lacy et al
#'
#' @param these_samples_metadata The metadata data frame that contains sample_id column with ids for the samples to be classified.
#' @param maf_data The MAF data frame to be used for matrix assembling. At least must contain the first 45 columns of standard MAF format.
#' @param seg_data The SEG data frame to be used for matrix assembling. Must be of standard SEG formatting, for example, as returned by get_sample_cn_segments.
#' @param sv_data The SV data frame to be used for matrix assembling. Must be of standard BEDPE formatting, for example, as returned by get_combined_sv.
#' @param projection The projection of the samples. Only used to retrerive data through GAMBLR when it is not provided. Defaults to grch37.
#' @param output The output to be returned after prediction is done. Can be one of predictions, matrix, or both. Defaults to both.
#' @param include_N1 Whether to set samples with NOTCH1 truncating mutations to N1 group as described in Runge et al (2021). Defaults to FALSE.
#' @return data frame with classification, binary matrix used in classification, or both
#' @import data.table randomForest dplyr readr stringr
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
        maf_data,
        recurrence_min = 3
    )

    # Mutations matrix
    lacy_feature_matrix$ssm <-
        get_coding_ssm_status(
            gene_symbols = lacy_features$ssm,
            these_samples_metadata = these_samples_metadata,
            maf_data = maf_data,
            include_hotspots = FALSE
        ) %>%
        column_to_rownames(
            "sample_id"
        )

    lacy_feature_matrix$ssm <- complete_missing_from_matrix(
        lacy_feature_matrix$ssm,
        these_samples_metadata$sample_id
    )

    # Hotspots matrix
    lacy_feature_matrix$hotspots <-
    maf_data %>%
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
        )

    colnames(lacy_feature_matrix$hotspots)[4:6] <-
        c(
            "POU2F2_239",
            "MYD88_265",
            "EZH2_646"
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

    lacy_feature_matrix$shm <-
    get_ashm_count_matrix(
        regions_bed = ashm_features,
        maf_data = maf_data,
        these_samples_metadata = these_samples_metadata
    )

    lacy_feature_matrix$shm <- complete_missing_from_matrix(
        lacy_feature_matrix$shm,
        these_samples_metadata$sample_id
    )

    lacy_feature_matrix$shm[lacy_feature_matrix$shm>0] = 1

    #  Amplifications were classed as driver events if they targeted a known
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
        as.data.table %>%
        setkey(
            chromosome,
            start,
            end
        ) %>%
        foverlaps(
            seg_data,
            .,
            nomatch = 0
        ) %>%
        dplyr::mutate(
            CN = 2*2^log.ratio,
            mutated = 1
        ) %>%
        # drop neutrals
        dplyr::filter(
            !data.table::between(
                CN,
                1,
                6,
                incbounds = FALSE
            )
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
        lacy_feature_matrix$ssm,
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
            str_detect(
              Variant_Classification,
              "Frame_Shift"
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
        )
      if (output == "both") {
      return(
        list(
          hmrn_matrix = lacy_feature_matrix$complete,
          hmrn_predictons = predictions
        )
      )
      }

    }

    if(output == "predictions"){
      return(predictions)
    }else if (output == "both") {
      return(
        list(
          lacy_matrix = lacy_feature_matrix$complete,
          lacy_predictons = predictions
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
