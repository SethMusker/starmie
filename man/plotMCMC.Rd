% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotMCMC.R
\name{plotMCMC}
\alias{plotMCMC}
\title{Plot STRUCTURE MCMC chains}
\usage{
plotMCMC(x, plot = TRUE, use_logL = TRUE, facet = TRUE)
}
\arguments{
\item{x}{\code{\link{structList}} objects or a \code{\link{struct}} object}

\item{plot}{logical print resulting plot default TRUE}

\item{use_logL}{logical plot log-likelihood (TRUE) or admixture coefficient}

\item{facet}{logical facet by K default TRUE}
}
\value{
If plot is TRUE a ggplot is printed to the screen and the
plot object and the data to generate it are returned. Otherwise,
a data.frame containing MCMC info it returned.
}
\description{
Plot non-burn MCMC iterations of STRUCTURE for checking convergence.
If plot is set to FALSE a data.frame is returned containing the log likelihood
and alpha values over different K and runs and not plot is printed to the device.
}
\examples{
#Read in Structure files
multiple_runs_k10 <- exampleStructure("mcmc_diagnostics")
print(multiple_runs_k10)
results <- plotMCMC(multiple_runs_k10, plot = TRUE)
single_run <- plotMCMC(multiple_runs_k10[[1]])
}
