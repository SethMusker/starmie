% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/averageQ.R
\name{averageQ}
\alias{averageQ}
\title{Convenience function for creating final Q-matrix after CLUMMP}
\usage{
averageQ(Q_list)
}
\arguments{
\item{Q_list}{A list of of Q matrices.}
}
\description{
Convenience function for creating final Q-matrix after CLUMMP
}
\examples{
# Read in Structure files
multiple_runs_k10 <- exampleStructure("mcmc_diagnostics")
Q_list <- lapply(multiple_runs_k10, getQ)
avgQ <- averageQ(Q_list)
}
