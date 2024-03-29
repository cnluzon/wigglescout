% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bwplot.R
\name{.violin_body}
\alias{.violin_body}
\title{Internal function to plot ranges in violin plot with a highlighted GRanges.}
\usage{
.violin_body(
  gr,
  highlight = NULL,
  minoverlap = 0L,
  highlight_label = NULL,
  highlight_colors = NULL
)
}
\arguments{
\item{gr}{GRanges object with as many columns as samples to plot.}

\item{highlight}{List of GRanges to use as highlight for subgroups.}

\item{minoverlap}{Minimum overlap required for a locus to be highlighted}

\item{highlight_label}{Labels for the highlight groups.
If not provided, column names are used.}

\item{highlight_colors}{Array of color values for the highlighting groups}
}
\value{
A named list where plot is a ggplot object and calculated is a list
of calculated values (for verbose mode).
}
\description{
Internal function to plot ranges in violin plot with a highlighted GRanges.
}
