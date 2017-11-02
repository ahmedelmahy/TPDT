\name{zhang_L2_boot}
\alias{zhang_L2_boot}
\title{L2-norm based test as described by Zhang in AoV of Functional Data - bootstrap version.}
\usage{
zhang_F(y0, t0, M = 100, alpha = 0.05, B = 1000, tmin, tmax)
}
\arguments{
\item{yo}{[list] containing all data of one trajectory for one individual per entry.}

\item{t0}{[list] containing the time points at which y0 were recorded with again one trajectory for one individual per entry..}

\item{M}{[Integer] number of points for which smooth functions are discretized.}

\item{alpha}{[numeric] significance level.}

\item{B}{[Integer] number of bootstrap realizations.}

\item{tmin}{[numeric] beginning of time interval for which differences should be evaluated with the test.}

\item{tmax}{[numeric] end of time interval for which differences should be evaluated with the test.}
}
\description{
Description to be added.
}

