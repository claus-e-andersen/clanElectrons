% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Bethe-stopping-power-computations-001.R
\name{electronic.MSP.restricted}
\alias{electronic.MSP.restricted}
\title{electronic.MSP.restricted}
\usage{
electronic.MSP.restricted(MeV = 1, delta.keV = 10, dat = NULL, delta = 0)
}
\arguments{
\item{MeV}{kinetic energy of the electron in MeV}

\item{delta.keV}{cut-off energy in keV}

\item{dat}{list with parameters}

\item{delta}{= density-effect correction}

\item{dat$Z}{= atomic number}

\item{dat$A}{= atomic mass}

\item{dat$I}{= mean excitation energy in eV}
}
\description{
Computation of restricted mass electronic stopping power (MSP)
according to ICRU-90 eq. 4.11 (page 22). The correction for the density
effect (delta) has to be supplied by the user in a seperate argument.
}
\details{
Notes:
The computed restricted MSP is in units of MeV pr. g/cm2.
}
