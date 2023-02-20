# clanElectrons (R package)
This package includes functions for computation of restricted and unrestricted mass electronic stopping power (Bethe formula) for electrons including Sternheimer density-effect corrections (both "exact" and by published parameters). The package is aimed for research applications within radiotherapy dosimetry and comparisons with ICRU-37 and ICRU-90 publications.

The density-effect correction accounts for the polarization of the material. Essentially, it explains why liquid water at energies above 550 keV or such has a mass electronic stopping power that is less than the mass electronic stopping power for water vapour at the same energy. The modelling involves dielectric stopping (i.e. optics or transport of electromagnetic waves in solid-state materials). The "stationary" electrons in the material are responsible for the absorption and re-emission of the "light" and the plasma frequency plays a vital role in the theory. The "stationary" electrons are modelled as oscillators with resonance frequencies determined by the electronic binding energies as given by the atomic configuration of the material. You may ask in what sense we can represent the electrical field from a single electron in some beam as an electromagnetic wave with a certain frequency? From the perspective of the material, a swift electron in a beam actually creates an electromagnetic field much like a delta-function in time. This field therefore can be represented as the superposition of plane waves with a very broad range of frequencies since the Fourier transformation of a delta function is unity. The density effect concerns only long-distance collisions (i.e. soft collisions). These gain more importance at high energies due to the Lorentz contraction of the field (it becomes pancake like at high energies rather than spherical).

Knowledge of the density effect can be of significant importance in composite materials like graphite or alanine dosimeters where the density controlling the "density effect" is the grain density rather than the bulk density of the material. Without proper accounting for these matters in Monte-Carlo modelling, one may wrongly introduce, for example, apparent changes in the intritic detector efficiency for high energies.   

# How to run it?
First, you will need R which can be downloaded from cran (https://cran.r-project.org/). Then you can download and run this R package (called clanElectrons) using the install_github("claus-e-andersen/clanElectrons") command from the devtools package. More details are given below. 

After a successful installation of the clanElectrons package, try to run the demonstration function:
```
demo.graphite.table() 
```
That should produce a table of computed stopping power data for graphite (discussed below).
To see the contents of the function (i.e. the function body), just call the function without "()":
```
demo.graphite.table
```
## Contact
https://orbit.dtu.dk/en/persons/claus-e-andersen

## Main references
- ICRU-90: Key data for ionizing-radiation dosimetry: Measurement standards and applications (2014/2016).
- R.M. Sternheimer (Brookhaven), M.J. Berger (NBS) and S.M. Seltzer (NBS): Density effect for the ionization loss of charged particles in various substances. Atomic Data and Nuclear Data Tables 30,26 l-27 1 (1984).
- G4DensityEffectCalculator.cc for Geant4 by Matthew Strait, School of Physics and Astronomy, University of Minneapolis, 2019. 
- T.A. Carlson: Photoelectron and Auger Spectroscopy (1975, Springer).

## Main functions in package
```
# The material in question is described in the list called 
#   dat 
# as will be explained below (see dat.H2O as an example).
# Some of the function calls (like Sternheimer.delta.exact) are designed
# to update dat with output such that dat both contains the material parameters
# and the result of the computation. 

# "Exact" Sternheimer density-effect correction (delta):
# See ICRU-90 equations 4.27, 4.28, and 4.29
dat <- Sternheimer.delta.exact(MeV = 1, dat)  

# Simple parametric model by Sternheimer for the density-effect correction (delta):
# See Sternheimer et al. 1984 paper
dat <- Sternheimer.delta.param(MeV = 1, dat) 


# Some functions in the package just output the results directly (like electronic.MSP).
# Here we do not update dat. 

# Electronic mass stopping power:
# See ICRU-90 equation 4.8
s <- electronic.MSP(MeV = 1, dat, delta = 0) 

# Electronic mass stopping power, restricted at delta.keV (i.e. the cut-off energy): 
# See ICRU-90 equation 4.11
s <- electronic.MSP.resticted(MeV = 1, delta.keV = 10, dat, delta = 0) 


# How to use Sternheimer.delta.exact and electronic.MSP:
demo.water(MeV = 1, delta.keV = 10)      
demo.graphite(MeV = 1, delta.keV = 10) 
demo.water.table()      
demo.graphite.table() 

# How to use Sternheimer.delta.param:
demo.Bragg.rule.test() 
```
The parameters (Z, A, I etc.) are in a list called dat, and the dentity-effect functions return a list of the same format. The output contains a copy of the parameters that were input to the function call plus whatever output was generated by the function. The MSP functions return a single number (i.e. the electronic mass stopping power).   

## Applications
You specify Z, A, I and other material parameters and the package has functions for computation
of electronic stopping power (both restricted and unrestricted) for electrons. The "exact" computation
of the density effect using Sternheimer theory requires knowledge of binding energies and occupation fractions for electrons in the
subshells of the material in question. 

We provide the material parameters to the functions using a list called dat (see dat.H2O and dat.graphite given below). 
Some of the important parameters for the "exact" Sternheimer model are dat$fvec, dat$Evec and dat$nc:

- fvec is a vector with electron occupation levels per subshell (number of electrons divided by the atomic number, Z). If the atom or coumpound is characterized by nlev sub-shells the this means that the material is modelled as having nlev oscillators.

- Evec is a vector of the same length as fvec with binding energies in eV. 

- nc is the number of conduction electrons per atom. For an insulator, nc = 0. In the example with graphite (see below), we set nc = 1. 

## Binding energies and sub-shell occupation levels
The atomic data required for the density effect computations can be found in the book by Carlson (1975, see above) starting page 337: Appendix 1
Atomic binding Energies for each subshell for elements  Z = 1-106.  ICRU-90 uses the Carlson data. The Geant4 implementation in G4DensityEffectCalculator.cc (see above) may have been based on other data.

## Solver details for the Sternheimer "exact" model
Two equations need to be solved in the function Sternheimer.delta.exact: one for the scaling parameter mu.st and one for L. A simple root-finding solution has been adapted where we map the functions over a relevant range of parameters and identify the regions where the functions contain both positive and negative values. We then find the root by linear regression.  The default parameters for the search algorithm are:

```
 mu.solver.parm <- list(
      mu.st.min = 0,
      mu.st.max = 20,
      mu.st.N   = 10000,
      mu.st.eps = 6e-05
    )

    L.solver.parm <- list(
      L.min = 0.02,
      L.max = 4000,
      L.N   = 80000,
      L.eps  = 1e-03
    )
```
where we provide minimum and maxium values to consider, the number of points where the function needs to be evaluated (N), and a criterium (eps) for how close to zero we need to be for the data included in the linear regression for finding the actual root. The user can supply alternative parameter sets in the call of Sternheimer.delta.exact. For example, we could limit the search for L to be between 0.02 and 20 as follows:
 ```
LSP <- list(
      L.min = 0.02,
      L.max = 20,
      L.N   = 80000,
      L.eps  = 1e-03
    )

dat <- Sternheimer.delta.exact(MeV = 1, dat, L.solver.parm = LSP) 
```
The default parameters should work in most cases, but tweeks may be needed in some situations. If roots are not found, that it may be useful to consult a plot of the functions. Such a plot can be retrieved from the function my setting: 
```
exact.plot = TRUE
```
in the dat-list. See the file plots001.pdf for an example of such a olot.


## Example 1: Stopping power computation for water and comparison with ICRU-90

First, we provide the data for water (liquid) in a list called dat.H2O:
```
dat.H2O <- list(
  Z    = 10,           # Atomic number
  A    = 18.0158,      # Atomic mass
  I    = 78,           # Mean excitation energy of medium in eV       
  exact.rho =  0.998,  # Density in g/cm3, needed for the "exact" density-effect correction
  nc   = 0,            # Number of conducting electrons pr. atom. Always treat compounds as insulators (i.e. nc = 0)
  fvec = c(2/10, 2/10, 2/10, 4/10),    # First 2 x H, then O
  Evec = c(13.6, 538.0, 28.48, 13.62), # Binding energies in eV for the subshells as given by Carlson (1975), see ICRU-90
  exact.plot = FALSE)
```

Secondly, we assign the parameter list to dat and do the computations:

```
dat <- dat.H2O

############################
# 800 keV
############################
MeV <- 0.8 # Electron kinetic energy

# Note that the results of the density-effect correction computation will be included in dat:
dat <- Sternheimer.delta.exact(MeV, dat)

# Compute MSP without density effect correction (delta = 0):
xx0 <- electronic.MSP.Bethe(MeV, dat, delta = 0) 

# Compute MSP with exact Sternheimer density correction:
xx <- electronic.MSP.Bethe(MeV, dat, delta = dat$exact.delta) 

# Collect the data in a dataframe. The following will produce the 800 keV row for the table shown below:
df1 <- data.frame(
MeV = MeV, 
I.eV = dat$I, 
rho = dat$exact.rho, 
MSP.R0 = xx0,
MSP.R = xx, 
MSP.ICRU90 = 1.880,
delta.R = dat$exact.delta, 
delta.ICRU90 = 0.1005,
mu.st = dat$mu.st,
L = dat$L)
```

## Excellent agreement between clanElectrons and ICRU-90 values for liquid water. 
```   
  MeV I.eV   rho     MSP.R0      MSP.R MSP.ICRU90     delta.R delta.ICRU90    mu.st            L
1e-03   78 0.998 118.059816 118.059816    118.100  0.00000000      0.00000 2.296822           NA
1e-02   78 0.998  22.384790  22.384790     22.390  0.00000000      0.00000 2.296822           NA
1e-01   78 0.998   4.092951   4.092951      4.093  0.00000000      0.00000 2.296822           NA
5e-01   78 0.998   2.025140   2.025140      2.025  0.00000000      0.00000 2.296822           NA
6e-01   78 0.998   1.958041   1.956419      1.956  0.01500487      0.01501 2.296822    0.5473502
8e-01   78 0.998   1.890531   1.880437      1.880  0.10044873      0.10050 2.296822    1.2677047
1e+00   78 0.998   1.864880   1.844806      1.845  0.20860797      0.20860 2.296822    1.7837140
1e+01   78 0.998   2.216852   1.966726      1.967  2.92799113      2.92800 2.296822   18.4868124
1e+02   78 0.998   2.798672   2.202281      2.202  6.99775891      6.99800 2.296822  195.1002875
1e+03   78 0.998   3.387159   2.400502      2.401 11.57725269     11.58000 2.296822 1957.7806954
```
where

  - MSP.R0 = mass electronic stopping power with delta = 0 (no correction of density effect).
  - MSP.R =  mass electronic stopping power computed with the clanElectrons software.
  - delta.R = the density-effect correction computed with the clanElectrons software.
  - Index ICRU90 = reference values from ICRU-90 (Appendix Table A.3).
  - mu.st = the scaling parameter from the oscillator strengths to I (ICRU-90 eq. 4.29).
  - L = the ell paramerer in ICRU-90 equation 4.28. 

The mass electronic stopping powers are given in units of MeV per g/cm2. 
See the functions demo.water() and demo.water.table() for further details.

Note that the density-effect correction (delta) for an insulator will be zero below a certain threshold.

The agreement is excellent in the sense that all four digits provided in ICRU-90 are replicated by the clanElectrons functions
if we allow for an uncertainty of +/- 1 on the last digit of the ICRU-90 results.

## Example 2: Stopping power computation for graphite and comparison with ICRU-90

First, we provide the data for graphite in a list called dat.graphite:
```
dat.graphite <- list(
    Z    = 6,       # Atomic number
    A    = 12.011,  # Atomic mass
    I    = 81,      # Mean excitation energy of medium in eV
    exact.rho =  2.265,          # Density in g/cm3, needed for the "exact" density-effect correction
    nc   = 1,                    # Number of conducting electrons pr. atom   
    fvec = c(2/6, 2/6, 1/6),     # Occupation fractions for the subshells in C
    Evec = c(288, 16.59, 11.26), # Binding energies in eV for the subshells as given by Carlson (1975), see ICRU-90
    exact.plot = FALSE           # Supplementary plots related to the root finding in the exact density correction
  )
```
These parameters are identical to what was used by Stefan Pojtinger and Ludwig Büermann from 
Physikalisch-Technische Bundesanstalt (PTB) in their paper "Characterization of new primary air kerma standards for dosimetry in Co-60, Cs-137 and Ir-192 gamma ray sources", Journal of Instrumentation, Volume 16, October 2021 (http://dx.doi.org/10.1088/1748-0221/16/10/P10014). Note that one electron per atom is modelled to be a "conducting, free electron" with zero binding energy whereas the remaining five electrons per atom are in three oscillators with given binding energies.  


Secondly, we assign the parameter list to dat and do the computations:

```
dat <- dat.graphite

############################
# 800 keV
############################
MeV <- 0.8 # Electron kinetic energy

# Note that the results of the density-effect correction computation will be included in dat:
dat <- Sternheimer.delta.exact(MeV, dat)

# Compute MSP without density effect correction (delta = 0):
xx0 <- electronic.MSP.Bethe(MeV, dat, delta = 0)               

# Compute MSP with exact Sternheimer density correction
xx  <- electronic.MSP.Bethe(MeV, dat, delta = dat$exact.delta) 

# Collect the data in a dataframe. The following will produce the 800 keV row for the table shown below:
  df1 <- data.frame(
    MeV = MeV,
    I.eV = dat$I,
    rho=dat$exact.rho,
    MSP.R0 = xx0,
    MSP.R = xx,
    MSP.ICRU90 = 1.640,
    delta.R = dat$exact.delta,
    delta.ICRU90 = 0.6075,
    mu.st = dat$mu.st,
    L = dat$L)
```


## Excellent agreement between clanElectrons and ICRU-90 values for graphite. 
```
    MeV I.eV   rho     MSP.R0      MSP.R MSP.ICRU90      delta.R delta.ICRU90    mu.st            L
  1e-03   81 2.265 104.766213 104.761358    104.800 2.470208e-04     0.000247 2.352228 2.557449e-02
  1e-02   81 2.265  19.993257  19.987943     19.990 2.634413e-03     0.002634 2.352228 8.185347e-02
  1e-01   81 2.265   3.664250   3.653922      3.654 4.047018e-02     0.040470 2.352228 2.937683e-01
  8e-01   81 2.265   1.694586   1.639646      1.640 6.074781e-01     0.607500 2.352228 1.693603e+00
  1e+00   81 2.265   1.671790   1.606030      1.606 7.593193e-01     0.759300 2.352228 2.069985e+00
  2e+00   81 2.265   1.694645   1.585502      1.586 1.364082e+00     1.364000 2.352228 3.834607e+00
  1e+01   81 2.265   1.989286   1.729347      1.729 3.381104e+00     3.381000 2.352228 1.838840e+01
  1e+02   81 2.265   2.512917   1.928154      1.928 7.623991e+00     7.624000 2.352228 1.962800e+02
  1e+03   81 2.265   3.042536   2.105602      2.106 1.221582e+01    12.220000 2.352228 1.957909e+03
```
where the symbols have the same meaning as in Example 1 except that now the ICRU-90 reference data are from
ICRU-90, Appendix Table A.2. For further details, see the functions demo.graphite() and demo.graphite.table().

Note that the density-effect correction (delta) for a conductor like graphite (where some electrons are essentially assigned zero
binding  energy), will have a non-zero value even for very low electron energies.

The agreement is excellent in the sense that all four digits provided in ICRU-90 are replicated by the clanElectrons functions
if we allow for an uncertainty of +/- 1 on the last digit of the ICRU-90 results.

## Installation in R or Rstudio

The library can be loaded into R using the install_github command which is in the devtools package. So you first need to ascertain that you have this package and you need to load it with the library command:

```
install.packages("devtools")
library(devtools)
install_github("claus-e-andersen/clanElectrons")
library(clanElectrons)

```
For the full functionalyty you will also need to load the packages lattice, dplyr and clanLattice:
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
```
You can also ask for help on any function, like so:
```
?clanElectrons
?demo.water
```
and then press on the index (at the bottom of the help page) to get the list of all functions. 
