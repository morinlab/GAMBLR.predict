% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classifiers.R
\name{complete_missing_from_matrix}
\alias{complete_missing_from_matrix}
\title{Complete samples missing from matrix.}
\usage{
complete_missing_from_matrix(
  incoming_matrix,
  list_of_samples,
  fill_in_values = 0,
  normalize_order = TRUE,
  samples_in_rows = TRUE
)
}
\arguments{
\item{incoming_matrix}{A matrix or data frame that should be filled. Required parameter.}

\item{list_of_samples}{Vector specifying all desired samples to be present in the resulting matrix. Required parameter.}

\item{fill_in_values}{Value that will be used to fill in the matrix.}

\item{normalize_order}{Logical parameter specifying whether sample order should be according to the supplied list. Default is TRUE.}

\item{samples_in_rows}{Logical argument indicating whether samples are in rows or columns. Default assumes samples are in rows and columns are features.}
}
\value{
A data frame with maintained orientation (rows and columns) where samples from the supplied list are present and reordered according to the specified order.
}
\description{
If some samples are missing from the matrix, add them with filled in 0 as value and normalize their ordering for consistency.
}
\examples{
partial_matrix = get_coding_ssm_status(these_samples_metadata = (get_gambl_metadata(case_set = "BL--DLBCL") \%>\% filter(pairing_status == "unmatched")), include_hotspots = FALSE)
complete_matrix = complete_missing_from_matrix(partial_matrix, get_gambl_metadata() \%>\% pull(sample_id))

}
