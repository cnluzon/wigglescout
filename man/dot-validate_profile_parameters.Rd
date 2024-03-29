% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{.validate_profile_parameters}
\alias{.validate_profile_parameters}
\title{Validate profile and heatmap relevant parameters}
\usage{
.validate_profile_parameters(bin_size, upstream, downstream)
}
\arguments{
\item{bin_size}{Bin size. Must be a positive number.}

\item{upstream}{Upstream bp. Must be positive and larger than bin size.}

\item{downstream}{Downstream bp. Must be positive and larger than bin size.}
}
\value{
Nothing (NULL) on success.
}
\description{
Validate profile and heatmap relevant parameters
}
