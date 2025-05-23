% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dlbclass.R
\name{construct_reduced_winning_version}
\alias{construct_reduced_winning_version}
\title{Construct reduced 21-dimension feature vector for DLBCLass}
\usage{
construct_reduced_winning_version(
  mutations = "inst/extdata/DLBCL.699.fullGSM.Sep_23_2022.tsv",
  fisher_test_result = "inst/extdata/fisher_exact_5x2.Sep_23_2022.combined.tsv",
  add_missing_features = FALSE
)
}
\arguments{
\item{fisher_test_result}{Fisher's test result from DLBCLass github repository}

\item{add_missing_features}{Set to TRUE to fill in missing features with zeroes}

\item{data}{Data frame of mutation status from DLBCLass supplement}
}
\value{
a list of data frames with the full mutation features, collapsed
21-dimension features and mutation-only (no CNV) features
}
\description{
This function is an R port of some of the pre-processing code
in DLBCLass https://github.com/getzlab/DLBCL-Classifier
}
\examples{

original_dlbclass_features = construct_reduced_winning_version()

}
