# clanElectrons (R package)
Electronic stopping power (Bethe formula) for electrons including Sternheimer density-effect corrections
(both "exact" and by parameters). Aimed for applications within radiotherapy dosimetry and comparisons with ICRU-37 
and ICRU-90 publications.

### Applications:
You specify Z, A, I, density and other material parameters and the package has functions for computation
of electronic stopping power (both respricted and unrestricted) for electrons. The "exact" computation
of the density effect using Sternheimer theory requires knowledge of binding energies for electrons in the
material in question.

### Sternheimer delta exact computation details:

We first supply material parameters:

  dat.Al.model1 <- list(
    plot.wanted = FALSE,
    MeV = 1000, # Kinetic energy
    nlev = 6,   # Number of subshells
    Z    = 13,  # Atomic number
    A    = 26.98154,      # Atomic mass
    rho.density =  2.265, # Density in g/cm3
    fvec.org = c(2/13, 2/13 ,2/13, 2/13, 2/13, 3/13), # Subshell occupancy level
    Evec.org = c(1564.0 , 121.0, 77.0, 77.0, 10.62, 5.986), # Binding energy for each subshell in eV
    I = 166.0 # Mean exicitation energy in eV
  )

Then we set the material to a conductor:
   
   dat.Al.model1 <- Sternheimer.set.to.conductor(dat.Al.model1)

or to an insulator:

   dat.Al.model2 <- Sternheimer.set.to.insulator(dat.Al.model1)

This involves manipulation of the binding energy of the outmost subshell (energy
will be set to zero for a conductor). 

Finally, compute the density correction factor (delta):

  dat.out1 <- Sternheimer.delta.exact(dat.Al.model1)
  dat.out2 <- Sternheimer.delta.exact(dat.Al.model2)

All parameters and output are kept in the dat.out lists.

## Installation in R or Rstudio

The library can be loaded into R using the install_github command which is in the devtools package:

install.packages("devtools")

library(devtools)

install_github("claus-e-andersen/clanElectrons")

library(clanElectrons)

To get a list of functions in the library, just call:

help(package=clanElectrons)

or

library(help=clanElectrons)
