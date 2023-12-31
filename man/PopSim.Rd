% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pop_sim.R
\name{PopSim}
\alias{PopSim}
\title{Virtual population simulation}
\usage{
PopSim(
  Amax = 10,
  Z = 0.45,
  k = 0.103,
  Linf = 63,
  t0 = -0.25,
  SD = 3,
  lmin = 0,
  lmax = 100,
  lbin = 0.5,
  Nyears = 20,
  Rmin = 1000,
  Rmax = 1e+06,
  Namos = 4
)
}
\arguments{
\item{Amax}{A value of the age max in the population}

\item{Z}{A value of total mortality (Z)}

\item{k}{A value of growth coeficient of the VBGF}

\item{Linf}{A value of asymptotic length}

\item{t0}{A value of t0 parameter of the VBGF}

\item{SD}{A value of the standard deviation of the length at age}

\item{lmin}{A value of minimum length of the P matrix}

\item{lmax}{A value of the maximum length of the P matrix}

\item{lbin}{A value of the length bins of the P matrix}

\item{Nyears}{Number of years that will be created}

\item{Rmin}{A value of the minimum Recruitment}

\item{Rmax}{A value of the maximum Recruitment}

\item{Namos}{Number of years that will be sampling}
}
\value{
Return a P matrix [rows = lengths, cols = ages].
}
\description{
It simulates a simple virtual population with specified characteristics
}
\examples{
Pkl(pars=c(60, 0.3), gpars = c(0, 3))

}
