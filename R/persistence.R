#' Load a previously saved DLBCLone model (including UMAP state)
#'
#' Reconstructs a DLBCLone model saved via [DLBCLone_save_optimized()] by
#' restoring the serialized model list and the associated `uwot` UMAP object.
#' Optionally performs an integrity check to verify that embeddings can be
#' reproduced from the saved state.
#'
#' @param path Character scalar. Directory containing the saved files.
#'   Defaults to \code{"."}.
#' @param name_prefix Character scalar. Common filename prefix used when the
#'   model was saved. The function looks for
#'   \code{<name_prefix>_model.rds} and \code{<name_prefix>_umap.uwot}
#'   in \code{path}. Defaults to \code{"DLBCLone"}.
#' @param check_integrity Logical. If \code{TRUE}, the function re-embeds the
#'   training samples (both batch and iterative modes, when available) and
#'   confirms that the coordinates match the saved embeddings (within a small
#'   tolerance). This requires that the saved model list contains
#'   \code{projection_train}, which is a consequence of activating the model
#'   with \code{DLBCLone_activate()}.
#' @param shiny_app_mode Logical. If \code{TRUE}, the function will just load
#'  as is and trust the existing pre-embeddings to speed up the load of shiny
#'  application.
#'
#' @return A list representing the reconstructed DLBCLone model, suitable for
#'   downstream use with \code{DLBCLone_predict()} or utilities that expect
#'   \code{DLBCLone_model}. On return, \code{DLBCLone_model$model} contains the
#'   restored \code{uwot} model object.
#'
#' @details
#' This function expects two files produced by \code{DLBCLone_save_optimized()}:
#' \itemize{
#' \item \code{<name_prefix>_model.rds}: a serialized list with the core model
#'       components (e.g., \code{features}, \code{df}, and optional test
#'       embeddings).
#' \item \code{<name_prefix>_umap.uwot}: a serialized \code{uwot} UMAP model
#'       saved with \code{uwot::save_uwot()}.
#' }
#'
#' If \code{check_integrity = TRUE}, the function computes fresh embeddings via
#' \code{make_and_annotate_umap()} and compares them to the saved ones. A
#' mismatch larger than \code{0.01} in either axis will throw an error.
#'
#' @seealso [DLBCLone_save_optimized()], \code{uwot::save_uwot()},
#'   \code{uwot::load_uwot()}
#'
#' @importFrom dplyr left_join mutate filter select
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Restore a model for use in a fresh R session
#' DLBCLone_myfeats_opt <- DLBCLone_load_optimized(
#'   path = ".",
#'   name_prefix = "DLBCLone_myfeats"
#' )
#'
#' # Restore and verify that saved embeddings can be reproduced
#' confirmed_DLBCLone <- DLBCLone_load_optimized(
#'   path = ".",
#'   name_prefix = "dlbclone_dlc_test",
#'   check_integrity = TRUE
#' )
#' }

DLBCLone_load_optimized <- function(
  path=".",
  name_prefix="DLBCLone",
  check_integrity = FALSE,
  shiny_app_mode = FALSE
){
  compare_embeddings = function(df1,df2){
    #match on sample_id
    df_comp = left_join(df1,df2,by="sample_id")
    df_comp = mutate(df_comp,delta_x = abs(V1.x - V1.y),
                   delta_y = abs(V2.x-V2.y))
    bad_rows = filter(df_comp,delta_x>0.01 | delta_y>0.01)
    if(nrow(bad_rows)){
      print(paste(nrow(bad_rows),"affected rows"))
      print(head(bad_rows))
      stop("Consistency issue detected!")
    }

  }
  prefix = paste0(path,"/",name_prefix)
  umap_file =  paste0(prefix,"_umap.uwot")

  rds_file=paste0(path,"/",name_prefix,"_model.rds")
  DLBCLone_model = readRDS(rds_file)
  umap_model = load_uwot(umap_file)
  DLBCLone_model$model = umap_model
  if(check_integrity){
    if("projection_train" %in% names(DLBCLone_model)){
      message("confirming integrity of model against existing embeddings using iterative mode")
      batch_test = make_and_annotate_umap(DLBCLone_model$features,
                                           DLBCLone_model$df %>% select(sample_id),
                                           DLBCLone_model,
                                           individually = TRUE)
      compare_embeddings(DLBCLone_model$projection_train,batch_test$df)
    }else{
      print(names(DLBCLone_model))
      stop("No training projection found in model. Cannot check integrity.")
    }
  }else{
    if (!shiny_app_mode){
        #remove potentially stale embeddings to force the user to re-project
        #if they want to use them. This ensures that only a model that has been verified
        #with check_integrity = TRUE will retain the embeddings.
        DLBCLone_model$projection_train = NULL
    }
  }
  return(DLBCLone_model)
}

#' Save a DLBCLone model (and optionally integrity test embeddings)
#'
#' Serializes the DLBCLone model list and its associated \code{uwot} UMAP
#' object to disk. Optionally computes and stores reference embeddings to enable
#' future integrity checks with [DLBCLone_load_optimized()].
#'
#' @param DLBCLone_model List. The model object to save. It must contain a
#'   trained UMAP model at \code{$model} created with \pkg{uwot}.
#' @param base_path Character scalar. Directory to write output files.
#'   Defaults to \code{"./"}.
#' @param name_prefix Character scalar. Prefix for the saved files. Files will
#'   be named \code{<name_prefix>_model.rds} and \code{<name_prefix>_umap.uwot}.
#'   Defaults to \code{"DLBCLone"}.
#' @param include_tests Logical. If \code{TRUE}, compute and embed two reference
#'   embeddings (batch and iterative modes) from the training data and store
#'   them in the saved model list as \code{embedding_batch} and
#'   \code{embedding_iterative}. Defaults to \code{FALSE}.
#' @param overwrite Logical. If \code{FALSE} (default), an error is thrown if
#'   either target file already exists. Set \code{TRUE} to replace existing
#'   files.
#'
#' @return Invisibly, the paths of the two written files (invisibly). Called for
#'   its side effect of writing \code{.rds} and \code{.uwot} files.
#'
#' @details
#' The UMAP model is stored separately using \code{uwot::save_uwot()} and then
#' removed from the list prior to writing the \code{.rds}, to avoid duplication.
#' When \code{include_tests = TRUE}, the function calls
#' \code{make_and_annotate_umap()} twice to create reproducible embeddings that
#' can be verified later with \code{check_integrity = TRUE} in
#' \code{DLBCLone_load_optimized()}.
#'
#' @seealso [DLBCLone_load_optimized()], \code{uwot::save_uwot()}
#'
#' @importFrom utils head
#' @export
#'
#' @examples
#' \dontrun{
#' # Save a trained model; error if files exist
#' DLBCLone_save_optimized(
#'   DLBCLone_model = trained_model,
#'   base_path = "save_optimized/trial_folder",
#'   name_prefix = "test_A"
#' )
#' }
#' \dontrun{
#' # Save with embedded test embeddings, overwriting if present
#' DLBCLone_save_optimized(
#'   DLBCLone_model = trained_model,
#'   base_path = "save_optimized/trial_folder",
#'   name_prefix = "test_A",
#'   include_tests = TRUE,
#'   overwrite = TRUE
#' )
#' }
DLBCLone_save_optimized = function(
    DLBCLone_model,
    base_path="./",
    name_prefix="DLBCLone",
    include_tests = FALSE,
    overwrite = FALSE
){

  umap_file =  paste0(base_path,"/",name_prefix,"_umap.uwot")
  rds_file=paste0(base_path,"/",name_prefix,"_model.rds")
  if(file.exists(umap_file) || file.exists(rds_file)){
    if(!overwrite){
      stop("The expected outputs already exist. Re-run with overwrite = TRUE to replace.")
    }else{
      message("Will replace existing files of the same name")
      unlink(umap_file)
      unlink(rds_file)
    }
  }
  if(include_tests){
    # ensure our model stores an example embedding with its training data so we can test that the same embedding can be recovered
    message("generating embeddings using batch mode")
    embedded_batch = make_and_annotate_umap(DLBCLone_model$features,
                                           DLBCLone_model$df %>% select(sample_id),
                                           DLBCLone_model,
                                           individually = FALSE)
    message("generating embeddings using iterative mode")
    embedded_iterative = make_and_annotate_umap(DLBCLone_model$features,
                                           DLBCLone_model$df %>% select(sample_id),
                                           DLBCLone_model,
                                           individually = TRUE)
    DLBCLone_model[["embedding_batch"]] = embedded_batch$df
    DLBCLone_model[["embedding_iterative"]] = embedded_iterative$df

  }
  # UMAP model must be saved separately using save_uwot
  umap_model = DLBCLone_model$model

  save_uwot(umap_model, umap_file)
  DLBCLone_model$model = NULL
  saveRDS(DLBCLone_model, rds_file)
  message("Saved model to ",rds_file," and UMAP model to ",umap_file)
}