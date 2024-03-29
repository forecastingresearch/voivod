% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metrics.R
\name{VoI_log}
\alias{VoI_log}
\title{Compute Log VoI, or expected KL divergence.}
\usage{
VoI_log(pu, puc, pc, punotc = NA, pc_weight = NA)
}
\arguments{
\item{pu:}{P(U)}

\item{puc:}{P(U|c)}

\item{pc:}{P(c)}

\item{punotc:}{P(U|¬c)}

\item{pc_weight:}{P(c) used for expected value calculation (defaults to
pc). Only use this if you want use a different P(c) for weighting than the
one provided by the forecaster.}
}
\description{
Compute Log VoI, or expected KL divergence.
}
