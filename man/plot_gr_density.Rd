% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bwplot.R
\name{plot_gr_density}
\alias{plot_gr_density}
\title{Density plot of a precalculated GRanges object}
\usage{
plot_gr_density(
  gr,
  x,
  y,
  plot_binwidth = 0.05,
  remove_top = 0,
  verbose = TRUE,
  selection = NULL
)
}
\arguments{
\item{gr}{Scored GRanges object}

\item{x}{Column in gr corresponding to the x axis}

\item{y}{Column in gr corresponding to the y axis.}

\item{plot_binwidth}{Resolution of the bins in the density histogram
(different to genomic bin size)}

\item{remove_top}{Return range 0-(1-remove_top). By default returns the
whole distribution (remove_top == 0).}

\item{verbose}{Put a caption with relevant parameters on the plot.}

\item{selection}{A GRanges object to restrict binning to a certain set of
intervals. This is useful for debugging and improving performance of
locus specific analyses.}
}
\value{
A ggplot object.
}
\description{
Plots a density plot from two mcols of a GRanges object files that were
already calculated using bw_loci or bw_bins
}
