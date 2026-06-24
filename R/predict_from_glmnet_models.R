#' @title Predict class labels using glmnet models from select_informative_features
#'
#' @description Apply fitted elastic-net logistic regression models (produced by
#'   \code{GAMBLR.utils::select_informative_features(return_models = TRUE)}) to
#'   a new feature matrix to predict class labels and/or class probabilities.
#'
#' @details
#' ## One-vs-rest models
#' Each model produces P(sample belongs to class X). Probabilities across
#' models are row-normalised so they sum to 1, and the class with the highest
#' normalised probability is the predicted label.
#'
#' ## Pairwise models
#' Each model predicts P(class A) over P(class B) for one pair. A
#' Copeland-style majority vote is used: for each sample, the class that
#' wins the most head-to-head comparisons is predicted. Ties are broken by
#' summing the win probabilities.
#'
#' ## Feature alignment
#' Each model was trained on a specific subset of features (stored in
#' \code{model$features}). Columns present in the model but absent from
#' \code{new_data} are filled with 0 (treated as wild-type). A warning is
#' issued if more than 10 \% of model features are missing.
#'
#' @param new_data A matrix or data frame with samples as rows and binary
#'   features as columns. Row names or a column named by \code{sample_id_column}
#'   are used as sample identifiers.
#' @param glmnet_output The list returned by
#'   \code{select_informative_features(return_models = TRUE)}. Must contain
#'   elements \code{models}, \code{contrast}, and \code{comparison_values}.
#' @param type Character. What to return:
#'   \describe{
#'     \item{\code{"class"}}{A data frame with columns \code{sample_id} and
#'       \code{predicted_class} (default).}
#'     \item{\code{"prob"}}{A data frame with \code{sample_id} and one column
#'       per class containing the class probability.}
#'     \item{\code{"all"}}{A data frame combining both \code{predicted_class}
#'       and per-class probability columns.}
#'   }
#' @param sample_id_column Column in \code{new_data} containing sample IDs.
#'   Set to \code{NULL} to use row names. Default \code{"Tumor_Sample_Barcode"}.
#'
#' @return A data frame. See \code{type} for column details.
#'
#' @import dplyr
#' @importFrom stats predict setNames
#' @export
#'
#' @examples
#' \dontrun{
#' mat  <- get_binary_matrix(these_samples_metadata = train_meta,
#'                           maf_data = train_maf)
#' out  <- GAMBLR.utils::select_informative_features(
#'   feature_matrix    = mat,
#'   metadata          = train_meta,
#'   comparison_column = "pathology",
#'   comparison_values = c("DLBCL", "FL", "BL"),
#'   return_models     = TRUE
#' )
#'
#' new_mat <- get_binary_matrix(these_samples_metadata = test_meta,
#'                              maf_data = test_maf)
#' preds <- predict_from_glmnet_models(new_mat, out)
#' preds_prob <- predict_from_glmnet_models(new_mat, out, type = "all")
#' }
predict_from_glmnet_models <- function(
    new_data,
    glmnet_output,
    type             = "class",
    sample_id_column = "Tumor_Sample_Barcode"
) {
    type <- match.arg(type, c("class", "prob", "all"))

    if (!requireNamespace("glmnet", quietly = TRUE)) {
        stop("Package 'glmnet' is required. Install with: install.packages('glmnet')")
    }

    # ── validate glmnet_output ────────────────────────────────────────────────
    required <- c("models", "contrast", "comparison_values")
    missing_fields <- setdiff(required, names(glmnet_output))
    if (length(missing_fields) > 0) {
        stop(
            "glmnet_output is missing elements: ",
            paste(missing_fields, collapse = ", "),
            ". Did you call select_informative_features(return_models = TRUE)?"
        )
    }

    models            <- glmnet_output$models
    contrast          <- glmnet_output$contrast
    comparison_values <- glmnet_output$comparison_values

    if (length(models) == 0) {
        stop("glmnet_output$models is empty; no models to predict with.")
    }

    # ── resolve sample IDs and feature matrix ─────────────────────────────────
    new_data <- as.data.frame(new_data)

    if (!is.null(sample_id_column) && sample_id_column %in% colnames(new_data)) {
        sample_ids <- new_data[[sample_id_column]]
        feat_cols  <- setdiff(colnames(new_data), sample_id_column)
        new_data   <- new_data[, feat_cols, drop = FALSE]
    } else {
        sample_ids <- rownames(new_data)
        if (is.null(sample_ids)) {
            stop(
                "new_data has no row names and sample_id_column '",
                sample_id_column, "' was not found. ",
                "Add row names or supply a valid sample_id_column."
            )
        }
    }

    # ── helper: build newx matrix for one model ───────────────────────────────
    make_newx <- function(model_entry) {
        needed  <- model_entry$features
        present <- intersect(needed, colnames(new_data))
        absent  <- setdiff(needed, colnames(new_data))

        if (length(absent) / length(needed) > 0.10) {
            warning(
                sprintf(
                    "Model '%s': %d of %d features absent from new_data (filled with 0). ",
                    model_entry$label, length(absent), length(needed)
                ),
                "Predictions may be unreliable."
            )
        } else if (length(absent) > 0) {
            message(
                sprintf(
                    "Model '%s': %d missing feature(s) filled with 0: %s",
                    model_entry$label, length(absent),
                    paste(absent, collapse = ", ")
                )
            )
        }

        mat <- matrix(0.0, nrow = nrow(new_data), ncol = length(needed),
                      dimnames = list(sample_ids, needed))
        if (length(present) > 0) {
            mat[, present] <- as.matrix(new_data[, present, drop = FALSE])
        }
        mat
    }

    # ── predict probabilities per model ───────────────────────────────────────
    n_samples <- nrow(new_data)

    if (contrast == "one_vs_rest") {
        # Each model: P(sample == pos_class)
        # comparison_values determines the column order in the prob matrix
        classes   <- comparison_values
        prob_mat  <- matrix(NA_real_, nrow = n_samples, ncol = length(classes),
                            dimnames = list(sample_ids, classes))

        for (cls in classes) {
            label <- paste(cls, "vs rest")
            mdl   <- models[[label]]
            if (is.null(mdl)) {
                warning("No model found for '", label, "'; class '", cls,
                        "' probabilities set to NA.")
                next
            }
            newx          <- make_newx(mdl)
            raw_prob      <- as.numeric(
                stats::predict(mdl$fit, newx = newx, s = mdl$lambda,
                               type = "response")
            )
            prob_mat[, cls] <- raw_prob
        }

        # Row-normalise so probabilities sum to 1 across classes
        row_sums <- rowSums(prob_mat, na.rm = TRUE)
        row_sums[row_sums == 0] <- 1  # avoid /0 for all-NA rows
        prob_norm <- sweep(prob_mat, 1, row_sums, "/")

        predicted_class <- classes[apply(prob_norm, 1, which.max)]
        prob_df <- as.data.frame(prob_norm)

    } else {
        # Pairwise: majority vote with probability tie-breaking
        classes  <- comparison_values
        # win_mat[i, cls] = total probability mass cls has "won" across all pairs
        win_mat  <- matrix(0.0, nrow = n_samples, ncol = length(classes),
                           dimnames = list(sample_ids, classes))

        for (mdl_name in names(models)) {
            mdl  <- models[[mdl_name]]
            if (is.null(mdl)) next
            newx <- make_newx(mdl)

            # P(pos_class wins this pair)
            p_pos <- as.numeric(
                stats::predict(mdl$fit, newx = newx, s = mdl$lambda,
                               type = "response")
            )
            p_neg <- 1 - p_pos

            pos_cls <- mdl$pos_class
            # The other class is derived from the label "A vs B"
            parts   <- strsplit(mdl$label, " vs ")[[1]]
            neg_cls <- parts[parts != pos_cls]
            if (length(neg_cls) != 1 || !neg_cls %in% classes) {
                warning("Cannot parse neg class from model label '", mdl$label,
                        "'; skipping.")
                next
            }

            if (pos_cls %in% classes) win_mat[, pos_cls] <- win_mat[, pos_cls] + p_pos
            if (neg_cls %in% classes) win_mat[, neg_cls] <- win_mat[, neg_cls] + p_neg
        }

        # Normalise win scores to probabilities
        row_sums <- rowSums(win_mat)
        row_sums[row_sums == 0] <- 1
        prob_norm       <- sweep(win_mat, 1, row_sums, "/")
        predicted_class <- classes[apply(prob_norm, 1, which.max)]
        prob_df         <- as.data.frame(prob_norm)
    }

    # ── assemble output ───────────────────────────────────────────────────────
    base <- data.frame(sample_id = sample_ids, stringsAsFactors = FALSE)

    if (type == "class") {
        base$predicted_class <- predicted_class
        return(base)
    }

    prob_df$sample_id <- sample_ids

    if (type == "prob") {
        return(dplyr::select(prob_df, sample_id, dplyr::everything()))
    }

    # type == "all"
    base$predicted_class <- predicted_class
    dplyr::left_join(base, prob_df, by = "sample_id")
}
