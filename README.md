# clanElectrons
Electronic stopping power (Bethe formula) for electrons including Sternheimer density-effect corrections
(both "exact" and by parameters). Aimed for comparison with ICRU-90.

The library can be loaded into R using the install_github command which is in the devtools package:

install.packages("devtools")

library(devtools)

install_github("claus-e-andersen/clanElectrons")

library(clanElectrons)

To get a list of functions in the library, just call:

help(package=clanElectrons)

or

library(help=clanElectrons)
