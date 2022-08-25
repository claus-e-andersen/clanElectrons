# clanElectrons (R package)
Restricted and unrestricted mass electronic stopping power (Bethe formula) for electrons including Sternheimer density-effect corrections
(both "exact" and by parameters). Aimed for research applications within radiotherapy dosimetry and comparisons with ICRU-37 
and ICRU-90 publications.

## Main references:
- ICRU-90: Key data for ionizing-radiation dosimetry: Measurement standards and applications (2014/2016).
- R.M. Sternheimer (Brookhaven), M.J. Berger (NBS) and S.M. Seltzer (NBS): Density effect for the ionization loss of charged particles in various substances. Atomic Data and Nuclear Data Tables 30,26 l-27 1 (1984).
- G4DensityEffectCalculator.cc for Geant4 by Matthew Strait (straitm-at-umn-dot-edu) 2019. 

## Applications:
You specify Z, A, I, density and other material parameters and the package has functions for computation
of electronic stopping power (both restricted and unrestricted) for electrons. The "exact" computation
of the density effect using Sternheimer theory requires knowledge of binding energies and occupation fractions for electrons in the
subshells of the material in question. fvec is a vector with electron occupation levels per subshell (number of electrons / Z).
Evec is a vector of the same length as fvec with binding energies in eV. nc is the number of conduction electrons per atom.
For an insulator, nc = 0. In the example with graphite (see below), we set nc = 1.

## Example 1: Stopping power computation for water and comparison with ICRU-90

First, we provide the data for water (liquid) in a list called dat.H2O:
```
dat.H2O <- list(
  Z    = 10,
  A    = 18.0158,
  I    = 78,
  exact.rho =  0.998,  # density in g/cm3
  nc   = 0,            # Number of conducting electrons pr. atom. Always treat compounds as insulators (i.e. nc =0)
  fvec = c(2/10, 2/10, 2/10, 4/10), # First 2 x H, then O
  Evec = c(13.6, 538.0, 28.48, 13.62),
  exact.plot = FALSE)
```

Notes:

exact referes to the detailed "exact" Sternheimer computation whereas param referes to the simplified 1984 model fits by Sternheimer. In the example below, we only use the exact method. We here include the param stuff for completeness. 

Secondly, we assign the parameter list to dat and du the computations:

```
dat <- dat.H2O

############################
# 800 keV
############################
MeV <- 0.8 # Electron kinetic energy
dat <- Sternheimer.delta.exact(MeV, dat)
xx0 <- electronic.MSP.Bethe(MeV, dat, delta = 0) # Compute MSP without density effect correction (delta = 0)
xx <- electronic.MSP.Bethe(MeV, dat, delta = dat$exact.delta) # Compute MSP with exact Sternheimer density correction

df1 <- data.frame(MeV = MeV, 
I.eV = dat$I, 
rho=dat$exact.rho, 
MSP.R0 = xx0,
MSP.R = xx, MSP.ICRU90 = 1.880,
delta.R = dat$exact.delta, 
delta.ICRU90 = 0.1005)
```

## Excellent agreement between clanElectrons and ICRU-90 values for liquid water. 
```
     MeV  I.eV   rho   MSP.R0    MSP.R MSP.ICRU90    delta.R delta.ICRU90
      0.8   78 0.998 1.890531 1.880437      1.880  0.1004488       0.1005
      1.0   78 0.998 1.864880 1.844806      1.845  0.2086075       0.2086
      10    78 0.998 2.216852 1.966726      1.967  2.9279906       2.9280
     100    78 0.998 2.798672 2.202281      2.202  6.9977588       6.9980
    1000    78 0.998 3.387159 2.400502      2.401 11.5772526      11.5800
```
where

  - MSP.R0 = mass electronic stopping power with delta = 0 (no correction of density effect).
  - MSP.R =  mass electronic stopping power computed with the clanElectrons software.
  - delta.R = the density-effect correction computed with the clanElectrons software.
  - Index ICRU90 = reference values from ICRU-90.

The mass electronic stopping powers are given in units of MeV per g/cm2 (i.e. approximately in MeV/cm since water has a density close to 1 g/cm3). 
See the function demo.Sternheimer.water() for further details.
  
## Example 2: Stopping power computation for graphite and comparison with ICRU-90

First, we provide the data for graphite in a list called dat.graphite:
```
dat.graphite <- list(
    Z    = 6,       # Atomic number
    A    = 12.011,  # Atomic mass
    I    = 81,      # Mean excitation energy in eV
    exact.rho =  2.265,         # Density in g/cm3, only needed for the exact density-effect correction.
    nc   = 1,                   # Number of conducting electrons pr. atom   
    fvec = c(2/6, 2/6, 1/6),    # Occupation fractions for the subshells in C
    Evec = c(288, 16.59,11.26), # Binding energies in eV of subshells from Carlson (1975), see ICRU-90.
    exact.plot = FALSE          # Supplementary plots related to the root finding in the exact density correction
  )
```
These parameters are identical to what was used by Stefan Pojtinger and Ludwig Büermann from 
Physikalisch-Technische Bundesanstalt (PTB) in their paper "Characterization of new primary air kerma standards for dosimetry in Co-60, Cs-137 and Ir-192 gamma ray sources", Journal of Instrumentation, Volume 16, October 2021 (http://dx.doi.org/10.1088/1748-0221/16/10/P10014). Note that one electron per atom is modelled to be a "conducting, free electron" with zero binding energy whereas the remaining five electrons per atom are in three oscillators with given binding energies.  


Secondly, we assign the parameter list to dat and du the computations:

```
dat <- dat.graphite

############################
# 800 keV
############################
MeV <- 0.8 # Electron kinetic energy
dat <- Sternheimer.delta.exact(MeV, dat)
xx0 <- electronic.MSP.Bethe(MeV, dat, delta = 0)               # Compute MSP without density effect correction (delta = 0)
xx  <- electronic.MSP.Bethe(MeV, dat, delta = dat$exact.delta) # Compute MSP with exact Sternheimer density correction

  df1 <- data.frame(
    MeV = MeV,
    I.eV = dat$I,
    rho=dat$exact.rho,
    MSP.R0 = xx0,
    MSP.R = xx,
    MSP.ICRU90 = 1.640,
    delta.R = dat$exact.delta,
    delta.ICRU90 = 0.6075)
```


## Excellent agreement between clanElectrons and ICRU-90 values for graphite. 
```
MeV   I.eV   rho   MSP.R0    MSP.R MSP.ICRU90    delta.R delta.ICRU90
  0.8   81 2.265 1.694586 1.639646      1.640  0.6074777       0.6075
  1.0   81 2.265 1.671790 1.606030      1.606  0.7593199       0.7593
  2.0   81 2.265 1.694645 1.585502      1.586  1.3640823       1.3640
  10    81 2.265 1.989286 1.729347      1.729  3.3811044       3.3810
 100    81 2.265 2.512917 1.928154      1.928  7.6239905       7.6240
1000    81 2.265 3.042536 2.105602      2.106 12.2158180      12.2200
```
where the symbols have the same meaning as in Example 1.

  
## Installation in R or Rstudio

The library can be loaded into R using the install_github command which is in the devtools package:

```
install.packages("devtools")
library(devtools)
install_github("claus-e-andersen/clanElectrons")
library(clanElectrons)

```
You will also need to load the packages lattice, dplyr and clanLattice:
```
library(lattice)
install.packages(dplyr)
library(dplyr)
install_github("claus-e-andersen/clanLattice")
library(clanLattice)
```

To get a list of functions in the library, just call one of these (after loading the package):
```
help(package=clanElectrons)
library(help=clanElectrons)
?clanElectrons
```
