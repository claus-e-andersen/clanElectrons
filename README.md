# clanElectrons (R package)
Electronic stopping power (Bethe formula) for electrons including Sternheimer density-effect corrections
(both "exact" and by parameters). Aimed for applications within radiotherapy dosimetry and comparisons with ICRU-37 
and ICRU-90 publications.

Applications:
You specify Z, A, I, density and other material parameters and the package has functions for computation
of electronic stopping power (both respricted and unrestricted) for electrons. The "exact" computation
of the density effect using Sternheimer theory requires knowledge of binding energies for electrons in the
material in question.

The library can be loaded into R using the install_github command which is in the devtools package:

install.packages("devtools")

library(devtools)

install_github("claus-e-andersen/clanElectrons")

library(clanElectrons)

To get a list of functions in the library, just call:

help(package=clanElectrons)

or

library(help=clanElectrons)
