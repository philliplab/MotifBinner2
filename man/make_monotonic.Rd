% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/consensus_cutoff.R
\name{make_monotonic}
\alias{make_monotonic}
\title{Makes a series of numbers monotonically increasing or decreasing}
\usage{
make_monotonic(x, direction = "increasing")
}
\arguments{
\item{x}{The vector of numbers to make monotonic}

\item{direction}{The direction in which the numbers should be made monotonic}
}
\description{
Monotonicity is achieved by substituting a number with the number preceeding
it if it violates the monotonicity criteria.
}
