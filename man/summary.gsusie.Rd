% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.gsusie.R
\name{summary.gsusie}
\alias{summary.gsusie}
\alias{print.summary.gsusie}
\title{Summarize G-SuSiE Fit.}
\usage{
\method{summary}{gsusie}(object, ...)

\method{print}{summary.gsusie}(x, ...)
}
\arguments{
\item{object}{A GSuSiE fit}

\item{\dots}{Additional arguments passed to the generic \code{summary}
or \code{print.summary} method.}

\item{x}{A (g)susie summary.}
}
\value{
a list containing a data frame of variables and a data frame
of credible sets.
}
\description{
\code{summary} method for the \dQuote{gsusie} class
Almost the same as in
\url{https://github.com/stephenslab/susieR/blob/master/R/summary.susie.R},
except that we change the output class into "summary.gsusie".
}
