#' load previously saved DLBCLone model and parameters
#' 
#' @param path Path to open saved the files
#' @param name_prefix Prefix of the saved files, all files will be in path and start with name_prefix
#'
#' @returns a reconstructed model for use with DLBCLone_predict
#' 
#' @import uwot
#' @import readr
#' @import dplyr
#' 
#' @export
#'
#' @examples
#' #Restore a model for use in a fresh R session
#' DLBCLone_myfeats_opt <- DLBCLone_load_optimized(
#'   path="./",
#'   name_prefix="DLBCLone_myfeats"
#' )
#'
DLBCLone_load_optimized <- function( 
  path=".",
  name_prefix="DLBCLone"
){
  prefix = paste0(path,"/",name_prefix)
  umap_file =  paste0(prefix,"_umap.uwot")

  rds_file=paste0(path,"/",name_prefix,"_model.rds")
  DLBCLone_model = readRDS(rds_file)
  umap_model = load_uwot(umap_file)
  DLBCLone_model$model = umap_model
  return(DLBCLone_model)
}


#' model storage for DLBCLone outputs
#' 
#' @param DLBCLone_model The model object to save
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
DLBCLone_save_optimized = function( 
    DLBCLone_model,
    base_path="./",
    name_prefix="DLBCLone"
){
  # UMAP model must be saved separately using save_uwot
  umap_model = DLBCLone_model$model
  umap_file =  paste0(base_path,"/",name_prefix,"_umap.uwot")
  save_uwot(umap_model, umap_file)
  DLBCLone_model$model = NULL
  rds_file=paste0(base_path,"/",name_prefix,"_model.rds")
  saveRDS(DLBCLone_model, rds_file)
message("Saved model to ",rds_file," and UMAP model to ",umap_file)
}