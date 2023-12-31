% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CPM_c.R
\name{CPM_c}
\alias{CPM_c}
\title{Convert FeatureCounts to CPM (counts per million reads)}
\usage{
CPM_c(counts_matrix, Assay)
}
\arguments{
\item{counts_matrix}{A Merged featureCounts table generated by \code{\link[=MergeFeatureCounts]{MergeFeatureCounts()}}.}

\item{Assay}{File prefix of the outputs".}
}
\description{
This function converts a Merged featureCounts table into CPM table. Merged FeatureCounts table is generated by
\code{\link[=MergeFeatureCounts]{MergeFeatureCounts()}}.
}
\details{
This function create two files: a table containing CPM in each sample
and a table containing average CPM among all samples.
}
