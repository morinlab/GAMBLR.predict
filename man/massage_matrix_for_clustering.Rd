% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classifiers.R
\name{massage_matrix_for_clustering}
\alias{massage_matrix_for_clustering}
\title{Will prepare the data frame of binary matrix to be used as NMF input. This means that for the features with SSM and CNV,
they will be squished together as one feature named GeneName-MUTorAMP or GeneName-MUTorLOSS, so the CNV features in the input data frame are expected
to be named GeneName_AMP or GeneName_LOSS. Next, for the genes with hotspot mutations labelled in the input data as
GeneNameHOTSPOT, the feature for hotspot mutation will be given preference and SSM with/without CNV will be set to 0 for that sample.
The naming scheme of the features as in this description is important, because the function uses regex to searh for these patters as specified.
Finally, if any features are provided to be dropped explicitly, they will be removed, and then the features not meeting the specified minimal
frequency will be removed, as well as any samples with 0 features.
Consistent with NMF input, in the input data frame each row is a feature, and each column is a sample. The input is expected to be numeric 1/0 with row and column names.}
\usage{
massage_matrix_for_clustering(
  incoming_data,
  blacklisted_cnv_regex = "3UTR|SV|HOTSPOT|TP53BP1|intronic",
  drop_these_features,
  min_feature_percent = 0.005
)
}
\arguments{
\item{incoming_data}{Input data frame or matrix to prepare for NMF.}

\item{blacklisted_cnv_regex}{Regular expression to match in feature names when considering SSM/CNV overlap.}

\item{drop_these_features}{Optional argument with features to drop from resulting matrix.}

\item{min_feature_percent}{Minimum frequency for the feature to be returned in the resulting matrix. By default, features present in less than 0.5\% of samples will be discarded.}
}
\value{
A matrix compatible with NMF input.
}
\description{
Will prepare the data frame of binary matrix to be used as NMF input. This means that for the features with SSM and CNV,
they will be squished together as one feature named GeneName-MUTorAMP or GeneName-MUTorLOSS, so the CNV features in the input data frame are expected
to be named GeneName_AMP or GeneName_LOSS. Next, for the genes with hotspot mutations labelled in the input data as
GeneNameHOTSPOT, the feature for hotspot mutation will be given preference and SSM with/without CNV will be set to 0 for that sample.
The naming scheme of the features as in this description is important, because the function uses regex to searh for these patters as specified.
Finally, if any features are provided to be dropped explicitly, they will be removed, and then the features not meeting the specified minimal
frequency will be removed, as well as any samples with 0 features.
Consistent with NMF input, in the input data frame each row is a feature, and each column is a sample. The input is expected to be numeric 1/0 with row and column names.
}
\examples{
data = system.file("extdata", "sample_matrix.tsv", package = "GAMBLR.predict") \%>\% read_tsv() \%>\% column_to_rownames("Feature")
NMF_input = massage_matrix_for_clustering(data)

}
