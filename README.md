# clanElectrons (R package)
Electronic stopping power (Bethe formula) for electrons including Sternheimer density-effect corrections
(both "exact" and by parameters). Aimed for applications within radiotherapy dosimetry and comparisons with ICRU-37 
and ICRU-90 publications.

## Applications:
You specify Z, A, I, density and other material parameters and the package has functions for computation
of electronic stopping power (both respricted and unrestricted) for electrons. The "exact" computation
of the density effect using Sternheimer theory requires knowledge of binding energies for electrons in the
material in question.

## Example: Stopping power computation for water and comparison with ICRU-90

First we provide the data for water (liquid):
```
dat.H2O <- list(
  Z    = 10,       # Atomic number
  A    = 18.0158,  # Atomic mass
  I    = 78,       # Mean excitation energy in eV
  #
  exact.rho =  0.998, # Density in g/cm3 needed for the exact density-correction
  exact.fvec = c(2/10, 2/10, 2/10, 4/10), # Occupation fractions for the subshells in H2 and O.
  exact.Evec = c(13.6, 538.0, 28.48, 13.62), # Binding energies of subshells from Carlson (1975), see ICRU-90.
  exact.plot = FALSE, # Supplementary plots related to the root finding in the exact density correction
  #
  param.note="Sternheimer et. al 1984, water (liquid) I = 75 and rho = 1.000 ", 
  param.C = -3.5017, param.X0 = 0.2400, param.X1 = 2.8004, param.a  = 0.09116, param.m  = 3.4773,
  param.delta.X0 = 0.097
  )
```

Notes:

exact referes to the detailed "exact" Sternheimer computation whereas param referes to the 1984 model fits by Sternheimer. In the example below we only use the exact method. We here include the param stuff for completeness. 

Always set compounds to insulators

```
dat.H2O <- Sternheimer.set.to.insulator(dat.H2O)
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

## Excellent agreement between the clanElectrons computations and ICRU-90 values for water. 

  MSP.R0 = mass electronic stopping power with delta = 0 (no correction of density effect).

  MSP.R =  mass electronic stopping power computed with the clanElectrons software.

  delta.R = the density-effect correction computed with the clanElectrons software.

  Index ICRU90 = reference values from ICRU-90.


     MeV I.eV   rho   MSP.R0    MSP.R MSP.ICRU90    delta.R delta.ICRU90
     0.80   78 0.998 1.890531 1.880437      1.880  0.1004488       0.1005
     1.00   78 0.998 1.864880 1.844806      1.845  0.2086075       0.2086
     10.0   78 0.998 2.216852 1.966726      1.967  2.9279906       2.9280
     100    78 0.998 2.798672 2.202281      2.202  6.9977588       6.9980
     1000   78 0.998 3.387159 2.400502      2.401 11.5772526      11.5800

## Installation in R or Rstudio

The library can be loaded into R using the install_github command which is in the devtools package:

```
install.packages("devtools")

library(devtools)

install_github("claus-e-andersen/clanElectrons")

library(clanElectrons)

```

To get a list of functions in the library, just call one of these (after loading the package):
```
help(package=clanElectrons)
library(help=clanElectrons)
?clanElectrons
```
