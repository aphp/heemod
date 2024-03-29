% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survival_operations.R
\name{apply_or}
\alias{apply_or}
\title{Apply an Odds Ratio}
\usage{
apply_or(dist, or, log_or = FALSE)
}
\arguments{
\item{dist}{A survival distribution.}

\item{or}{An odds ratio to be applied.}

\item{log_or}{If \code{TRUE}, the odds ratio is exponentiated
before being applied.}
}
\value{
A \code{surv_po} object.
}
\description{
Proportionally increase or reduce the odds of an event of
a survival distribution.
}
\examples{

dist1 <- define_surv_dist(distribution = "exp", rate = .25)
po_dist <- apply_or(dist1, 1.2)
}
