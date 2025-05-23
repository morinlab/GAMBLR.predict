% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/umap.R
\name{DLBCLone_train_test_plot}
\alias{DLBCLone_train_test_plot}
\title{Plot the result of a lymphgen classification}
\usage{
DLBCLone_train_test_plot(
  test_df,
  train_df,
  predictions_df,
  other_df,
  details,
  annotate_accuracy = FALSE,
  classes = c("BN2", "ST2", "MCD", "EZB", "N1"),
  label_offset = 2,
  title1 = "GAMBL",
  title2 = "predicted_class_for_HighConf",
  title3 = "predicted_class_for_Other"
)
}
\arguments{
\item{test_df}{Data frame containing the test data with UMAP coordinates}

\item{train_df}{Data frame containing the training data with UMAP coordinates}

\item{predictions_df}{Data frame containing the predictions with UMAP coordinates}

\item{other_df}{Data frame containing the predictions for samples in the "Other" class}

\item{details}{Single-row data frame with the best parameters from DLBCLone_optimize_params}

\item{annotate_accuracy}{Set to true to add labels with accuracy values}

\item{classes}{Vector of classes that were used in the training and testing}

\item{label_offset}{Length of the label offset for the accuracy labels}

\item{title1}{additional argument}

\item{title2}{additional argument}

\item{title3}{additional argument}
}
\value{
a ggplot object
}
\description{
Plot the result of a lymphgen classification
}
\examples{
#add the dataset name to the metadata if it's not already there (required for the plot work)
lymphgen_A53_DLBCLone$df$dataset = "GAMBL"

DLBCLone_train_test_plot(
 test_df = lymphgen_A53_DLBCLone$df,
 train_df = lymphgen_A53_DLBCLone$df,
 predictions_df = lymphgen_A53_DLBCLone$predictions,
 #other_df = lymphgen_A53_DLBCLone$predictions_other, #required only when "Other" was in the truth_classes
 details = lymphgen_A53_DLBCLone$best_params,
 classes = c("MCD","EZB","BN2","ST2","N1","A53","Other"),
 annotate_accuracy=TRUE,label_offset = 1)

}
