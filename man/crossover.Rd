% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Meiosis.R, R/RcppExports.R
\name{crossover}
\alias{crossover}
\alias{crossover}
\title{Simulate crossover locations.}
\usage{
crossover(L, m, p, obligate_chiasma, Lstar)
}
\arguments{
\item{L}{Double. Length of the chromosome in cM.}

\item{m}{Integer. The interference parameter (\code{m = 0} is no interference).}

\item{p}{Double. Proportion of chiasmata from no-interference mechanism.
(\code{p = 0} gives pure chi-square model)}

\item{obligate_chiasma}{Logical. If TRUE, require an obligate chiasma on the
4-strand bundle at meiosis. Only possible if all chromosomes are longer than 50 cM.}

\item{Lstar}{Double. Reduced chromosome length as produced by \code{\link{calc_Lstar}}.}
}
\value{
Double Vector. Crossover locations.
}
\description{
Simulate crossover locations on a single meiotic product using the
Stahl model.
}
\details{
This function is an R-wrapper of an underlying C++ routine.
It is not intended for direct usage, but exposed for completeness.
}
\section{Model}{

Chiasma locations are a superposition of two
processes: a proportion p exhibiting no interference, and a
proportion \code{(1 - p)} following the chi-square model with interference
parameter m.  Crossover locations are derived by thinning the
chiasma locations with probability \code{1/2}.

Simulations are under the Stahl model with the
interference parameter being an integer. This is an extension of
the chi-square model, but with chiasmata being the superposition of
two processes, one following the chi-square model and the other
exhibiting no interference.
}

\examples{
Meiosis::crossover(300, 10, 0.5, FALSE, 300)

}
\references{
Copenhaver, G. P., Housworth, E. A. and Stahl, F. W. (2002) Crossover
interference in arabidopsis.  \emph{Genetics} \bold{160}, 1631--1639.

Foss, E., Lande, R., Stahl, F. W. and Steinberg, C. M. (1993) Chiasma
interference as a function of genetic distance. \emph{Genetics}
\bold{133}, 681--691.

Zhao, H., Speed, T. P. and McPeek, M. S. (1995) Statistical analysis
of crossover interference using the chi-square model.  \emph{Genetics}
\bold{139}, 1045--1056.
}
