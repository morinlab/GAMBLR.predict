#' BL Classifier model.
#'
#' A random forest model developed and described in PMID: 36201743. Used to
#' classify BL tumours into genetic subgroups.
#'
#' @format ## `RFmodel_BL`
#' A randomForest model.
"RFmodel_BL"

#' FL Classifier model.
#'
#' A random forest model developed and described in PMID: 37084389. Used to
#' classify FL tumours into genetic subgroups.
#'
#' @format ## `RFmodel_FL`
#' A randomForest model.
"RFmodel_FL"

#' DLBCL Classifier model.
#'
#' A random forest model developed to reproduce the DLBCL grouping approach
#' described in PMID: 32187361. Used to classify DLBCL tumours into HMRN groups.
#'
#' @format ## `RFmodel_Lacy`
#' A randomForest model.
"RFmodel_Lacy"

#' Features for DLBCL grouping by Chapuy method.
#'
#' A named list of data frames and lists needed to assemble matrix and classify
#' DLBCL tumors according to Chapuy method described in PMID: 29713087.
#'
#' @format ## `chapuy_features`
#' A named list.
"chapuy_features"

#' Features for DLBCL grouping by Lacy method.
#'
#' A named list of data frames and lists needed to assemble matrix and classify
#' DLBCL tumors according to Lacy method described in PMID: 32187361.
#'
#' @format ## `lacy_features`
#' A named list.
"lacy_features"

#' Features for DLBCL grouping by unified method.
#'
#' A named list of data frames and lists needed to assemble matrix and classify
#' DLBCL tumors according to the unified classification.
#'
#' @format ## `lymphgenerator_features`
#' A named list.
"lymphgenerator_features"

#' All DLBCLass features.
#'
#' The features associated with each DLBCLass class according to the Supplemental
#'      data from Chapuy et al, 2024.
#'
#' @format ## `dlbclass_features_all`
#' A data frame with 163 rows and 7 columns.
#' \describe{
#'   \item{gene}{Gene symbol in Hugo format}
#'   \item{cluster}{The associated cluster}
#'   \item{AlterationType}{The type of alteration}
#'   \item{Marker.Gene.in.Cluster}{TODO}
#'   \item{Unique}{TODO}
#'   \item{Rank.within.Cluster}{Rank of feature importance for the class}
#'   \item{q}{q-value}
#' }
#' @keywords internal
"dlbclass_features_all"
