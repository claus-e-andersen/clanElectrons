# Created: July 28, 2022
# Revised: July 31, 2022
# Revised: August 28, 2022
# Name   : Claus E. Andersen

# Mass electronic stopping Power for electrons
# Restricted mass electronic stopping Power for electrons
# Sternheimer "exact" modelling of the density-effect correction
# Use of Sternheimer's fitter parameters for the density-effect correction

#' @title electronic.MSP
#' @description  Computation of mass electronic stopping power (MSP) according to
#' ICRU-90 eq. 4.8 (page 22). The correction for the density effect (delta)
#' has to be supplied by the user in a separate argument.
#'
#' dat$Z = atomic number
#' dat$A = atomic mass
#' dat$I = mean excitation energy in eV
#' delta = density-effect correction
#'
#' @param MeV kinetic energy of the electron in MeV
#' @param dat list with parameters
#' @param delta delta parameter
#' @return mass stopping power
#' @details Notes:
#' The computed MSP is in units of MeV pr. g/cm2.
#' @export
############################################################################
# Mass electronic stopping power
############################################################################
electronic.MSP <- function(MeV = 1, dat = NULL(), delta=0){
# Created: July 29, 2022
# Revised: August 27, 2022
# Name:    Claus E. Andersen
# Input:
#   MeV = kinetic energy of electron in MeV (this can be a vector)
#   I   = ionization energy in eV
#   Z   = atomic number (e.g. 6 for graphite)
#   A   = atomic mass (e.g. 12.011 for graphite)
# Output:
#   The mass electronic stopping power for electrons in units of MeV pr. g/cm2

T <- MeV
Z <- dat$Z
A <- dat$A
I <- dat$I

re     <- 2.81794092e-15 # m, classical electron radius
e      <- 1.60217733e-19 # J
E0.MeV <- 0.51099895000
E0     <- E0.MeV  * 1E6 * e # J
Na     <- 6.0221367e+23
u      <- 1.6605402e-27 # kg, atomic mass unit = 1/Na * 1000
c      <-  299792458 # m/s
tau    <- T / E0.MeV

beta <- (1 - (E0.MeV/(E0.MeV+T))^2 )^0.5
Fminus <- (1-beta^2)*(1 + tau^2/8 - (2*tau+1)*log(2))

P1 <- 2*pi*re^2*E0/beta^2 * Z * Na /  A
P2 <- 2*log(T*1e6/I) + log(1+tau/2) + Fminus - delta
MSP <- P1 * P2 /( 1E6 * e) * 10000 # MeV pr. g/cm2

MSP
} # End function


#' @title electronic.MSP.restricted
#' @description  Computation of restricted mass electronic stopping power (MSP)
#' according to ICRU-90 eq. 4.11 (page 22). The correction for the density
#' effect (delta) has to be supplied by the user in a seperate argument.
#'
#' dat$Z = atomic number
#' dat$A = atomic mass
#' dat$I = mean excitation energy in eV
#' #'
#' @param MeV kinetic energy of the electron in MeV
#' @param delta.keV cut-off energy in keV
#' @param dat list with parameters
#' @param delta = density-effect correction
#' @return restricted electronic mass stopping power

#' @details Notes:
#' The computed restricted MSP is in units of MeV pr. g/cm2.
#' @export

############################################################################
# Restricted Mass electronic stopping power
############################################################################
electronic.MSP.restricted <- function(MeV = 1, delta.keV = 10, dat=NULL, delta=0){
# Created: July 29, 2022
# Revised: July 31, 2022
# Revised: Aug 27, 2022
# Name:    Claus E. Andersen
# Input:
#   MeV = kinetic energy of electron in MeV (this can be a vector)
#   I   = ionization energy in eV
#   delta.kV = threshold energy in keV for the restriction
#   Z   = atomic number (e.g. 6 for graphite)
#   A   = atomic mass (e.g. 12.011 for graphite)
# Output:
#   The restricted mass electronic stopping power for electrons in units of MeV pr. g/cm2

# The max delta is 50% of the kinetic energy
delta.keV <- pmin(delta.keV, 0.5*MeV*1000)

T      <- MeV
Z      <- dat$Z
A      <- dat$A
I      <- dat$I

re     <- 2.81794092e-15 # m, classical electron radius
e      <- 1.60217733e-19 # J
E0.MeV <- 0.51099895000
E0     <- E0.MeV  * 1E6 * e # J
Na     <- 6.0221367e+23
u      <- 1.6605402e-27 # kg, atomic mass unit = 1/Na * 1000
c      <-  299792458 # m/s
tau    <- T / E0.MeV
eta    <- delta.keV / 1000 / T
beta   <- (1 - (E0.MeV/(E0.MeV+T))^2 )^0.5

# ICRU-90 eq. 4.12
Hminus <-  -1 - beta^2 + log( 4*(1-eta)*eta) + (1-eta)^-1 + (1-beta^2)*(tau^2*eta^2/2 + (2*tau+1) * log(1-eta))

P1 <- 2*pi*re^2*E0/beta^2 * Z * Na /  A
P2 <- 2*log(T*1e6/I) + log(1+tau/2) + Hminus - delta
rMSP <- P1 * P2 /( 1E6 * e) * 10000 # MeV pr. g/cm2
rMSP
} # End function


#' @title Sternheimer.data.read
#' @description  Data from Sternheimers approximitive fit to the density-effect correction
#' in the paper:
#' R.M. Sternheimer (Brookhaven), M.J. Berger (NBS) and S.M. Seltzer (NBS): Density effect
#' for the ionization loss of charged particles in various substances. Atomic Data and Nuclear
#' Data Tables 30,26 l-27 1 (1984).
#' @param print.wanted = TRUE or FALSE
#
#' @details
#' The function returns a data frame with fitted parameters for a few selected materials.
#' @export

############################################################################
# Sternheimer model for density effect: Data for different materials
############################################################################
Sternheimer.data.read <- function(print.wanted=FALSE){
# Created: July 29, 2022
# Revised: July 31, 2022
# Name:    Claus E. Andersen
# Data from the paper:
#   DENSITY EFFECT FOR THE IONIZATION LOSS OF CHARGED PARTICLES
#   IN VARIOUS SUBSTANCES
#   by R. M. STERNHEIMER (Brookhaven), M. J. BERGER (NBS) and S. M. SELTZER (NBS).
#   ATOMIC DATA AND NUCLEAR DATA TABLES 30,26 l-27 1 ( 1984)
#
df <- rbind(
data.frame(
Z=1,
Sternheimer.tab.id = "hydrogen-1",
note="Sternheimer et. al 1984, I = 19.2 rho = 8.3748E-5",
C = -9.5835, X0 = 1.8639, X1 = 3.2718, a  = 0.14092, m  = 5.7273, delta.X0 = 0.00
),

data.frame(
Z=1,
Sternheimer.tab.id = "hydrogen-2",
note="Sternheimer et. al 1984, I = 21.8 rho = 6.0000e-2 (liquid)",
C = -3.2632, X0 = 0.4759, X1 = 1.9215, a  = 0.13483, m  = 5.649, delta.X0 = 0.00
),

data.frame(
Z=6,
Sternheimer.tab.id = "graphite-1",
note="Sternheimer et. al 1984, I = 78 and rho = 2.265",
C = -2.8680, X0 = -0.0178, X1 = 2.3415, a  = 0.26142, m  = 2.8697, delta.X0 = 0.12
),

data.frame(
Z=6,
Sternheimer.tab.id = "graphite-2",
note="Sternheimer et. al 1984, I = 78 and rho = 2",
C = -2.9925, X0 = -0.0351, X1 = 2.486, a  = 0.20240, m  = 3.0036, delta.X0 = 0.10
),

data.frame(
Z=6,
Sternheimer.tab.id = "graphite-3",
note="Sternheimer et. al 1984, I = 78 and rho = 1.7",
C = -3.1550, X0 = 0.0480, X1 = 2.5387, a  = 0.20762, m  = 2.9532, delta.X0 = 0.14
),

data.frame(
Z=8,
Sternheimer.tab.id = "oxygen",
note="Sternheimer et. al 1984, I = 95 and rho = 1.3315e-3",
C = -10.7004, X0 = 1.7541, X1 = 4.3213, a  = 0.11778, m  = 3.2913, delta.X0 = 0.00
),

data.frame(
Z=82,
Sternheimer.tab.id = "lead",
note="Sternheimer et. al 1984, I = 823 and rho = 11.35",
C = -6.2018, X0 = 0.3776, X1 = 3.8073, a  = 0.09359, m  = 3.1608, delta.X0 = 0.14
),

data.frame(
Z=10,
Sternheimer.tab.id = "water-1",
note="Sternheimer et. al 1984, I = 75 and rho = 1.000 (liquid)",
C = -3.5017, X0 = 0.2400, X1 = 2.8004, a  = 0.09116, m  = 3.4773, delta.X0 = 0.097
),

data.frame(
Z=10,
Sternheimer.tab.id = "water-2",
note="Sternheimer et. al 1984, I = 71.6 and rho = 7.5618E-4 (vapor)",
C = -10.5962, X0 = 1.7952, X1 = 4.3437, a  = 0.08101, m  = 3.5901, delta.X0 = 0.121
)

)
if(print.wanted){
print("Reading Sternheimer data")
print(df)
}
df
}# end function


#' @title Sternheimer.delta.param
#' @description  Compute the density-effect correction based on the fitting parameters
#' @details
#' The fitted parameters can be found, for example, in the paper:

#' - R.M. Sternheimer (Brookhaven), M.J. Berger (NBS) and S.M. Seltzer (NBS): Density effect
#' for the ionization loss of charged particles in various substances. Atomic Data and Nuclear
#' Data Tables 30,26 l-27 1 (1984).
#'
#' The function returns an approximate value for the density-effect correction (delta).
#' See also Andero et al. "Fundamentals of ionizing radiation dosimetry" (2017) p. 74.
#'
#' The output ia returned as a list (dat):
#'   dat$param.MeV =  MeV
#'   dat$param.delta = delta
#'
#' First create a material (i.e. make a dat-list) called dat or something like than.
#'
#' Then call the function for the given energy MeV:
#'
#' dat <- Sternheimer.delta.param(MeV, dat)
#'
#' dat$param.C = Sternheimer parameter (see 1984 paper)
#' dat$param.X0 = Sternheimer parameter (see 1984 paper)
#' dat$param.X1 = Sternheimer parameter (see 1984 paper)
#' dat$param.a = Sternheimer parameter (see 1984 paper)
#' dat$param.m = Sternheimer parameter (see 1984 paper)

#'
#' @param MeV = kinetic energy of the electron (in MeV)
#' @param dat = list with parameters
#' @return a list (dat)
#' @export
############################################################################
# Sternheimer model for density effect
############################################################################
Sternheimer.delta.param <- function(MeV = 1, dat = NA){
# Created: July 29, 2022
# Revised: August 26, 2022
# Name:    Claus E. Andersen
# Input:
#   MeV = kinetic energy of electron in MeV (this can be a vector)
#   I   = ionization energy in eV
#   Z   = atomic number (e.g. 6 for graphite)
#   A   = atomic mass (e.g. 12.011 for graphite)
#   param.C, paran.X0 etc. : Sternheimer model fit parameters
#   # Output:
#   The Sternheimer density effect correction (delta) needed in the Bethe formula.
E0   <- 0.51099895000
beta <- (1 - (E0/(E0+MeV))^2 )^0.5
X  <- log10(beta / (1-beta^2)^0.5)

# Sternheimer parameters
C  <- dat$param.C
X0 <- dat$param.X0
X1 <- dat$param.X1
a  <- dat$param.a
m  <- dat$param.m
delta.X0 <- dat$param.delta.X0

# Sternheimer models for different energy regimes:
delta0 <- 10^(2*(X-X0))*delta.X0
delta1 <- 4.6052 * X + a*(X1-X)^m + C
delta2 <- 4.6052 * X  + C

delta <- delta0

ok <- (X0 < X) & (X < X1)
if(sum(ok)>0){
  delta[ok] <- delta1[ok]
}

ok <- X > X1
if(sum(ok)>0){
  delta[ok] <- delta2[ok]
}

dat$param.MeV <- MeV
dat$param.delta <- delta

# Return the density effect correction, delta.
dat
} # Sternheimer.delta.param

#' @title demo.Sternheimer.delta.param
#' @description  Demonstration of how to use Sternheimer.delta.param().
#' @details
#' None
#' @export
############################################################################
# Stopping power validation computations
############################################################################
demo.Sternheimer.delta.param <- function(){
# Created: July 29, 2022
# Revised: August 26, 2022
# Name:    Claus E. Andersen
print("This function should be run manually, line by line.")

dat.C <- list(
    Z    = 6,
    A    = 12.011,
    I    = 81,
    param.note="Sternheimer et. al 1984, graphite w. I = 78 and rho = 2.265",
    param.C = -2.8680, param.X0 = -0.0178, param.X1 = 2.3415, param.a  = 0.26142, param.m  = 2.8697,
    param.delta.X0 = 0.12
    )
dat <- dat.C

# Note that the param data are for I = 78  whereas the Bethe eq. is evaluated at I = 81
# Sternheimer parameter model
MeV <- 0.01
out <- Sternheimer.delta.param(MeV,dat)
df1 <- data.frame(MeV=out$param.MeV, delta=out$param.delta, delta.ICRU=0.002634)

MeV <- 0.1
out <- Sternheimer.delta.param(MeV,dat)
df2 <- data.frame(MeV=out$param.MeV, delta=out$param.delta, delta.ICRU=0.04047)

MeV <- 1
out <- Sternheimer.delta.param(MeV,dat)
df3 <- data.frame(MeV=out$param.MeV, delta=out$param.delta, delta.ICRU=0.7593)

MeV <- 10
out <- Sternheimer.delta.param(MeV,dat)
df4 <- data.frame(MeV=out$param.MeV, delta=out$param.delta, delta.ICRU=3.381)

rbind(df1,df2,df3,df4)
} # demo.Sternheimer.delta.param

#' @title demo.electronic.MSP.restricted
#' @description  Demonstration of how to use some the restricted electronic stopping power.
#' @details
#' Here we look at water with cutoff values equal to 1, 10, 100, and 1000 keV.
#' @export
############################################################################
# Stopping power demonstration
############################################################################
demo.electronic.MSP.restricted <- function(){
# Created: July 29, 2022
# Revised: July 31, 2022
# Revised: Aug 27, 2022
# Name:    Claus E. Andersen
print("This function should be run manually, line by line.")

    dat.H2O <-
    list(
      I=75,
      Z=10,
      A=18.0158)

    dat <- dat.H2O

xx <- seq(-3,2,length=50)
ee <- 10^xx
yy1 <- electronic.MSP.restricted(ee, delta.keV=1, dat, delta=0)
yy2 <- electronic.MSP.restricted(ee, delta.keV=10, dat, delta=0)
yy3 <- electronic.MSP.restricted(ee, delta.keV=100, dat, delta=0)
yy4 <- electronic.MSP.restricted(ee, delta.keV=1000, dat, delta=0)

df <- rbind(
      data.frame(MeV=ee, Sel.rho = yy1, keV.delta=1),
      data.frame(MeV=ee, Sel.rho = yy2, keV.delta=10),
      data.frame(MeV=ee, Sel.rho = yy3, keV.delta=100),
      data.frame(MeV=ee, Sel.rho = yy4, keV.delta=1000)
)

lattice::xyplot(log10(Sel.rho) ~ log10(MeV),
main="Restricted mass electronic stopping power for water (I = 75 eV)",
auto.key=list(title="Delta [keV]",columns=4),
groups=keV.delta,
data=df)
}


#' @title demo.Bragg.rule.test
#' @description  This function demonstrates the application of Bragg's rule
#' @details
#' We compute the stopping power for water based on the stopping powers of H and O and we
#' then compare that with a computation for H2O in one go.
#'
#' What we can learn from this function:
#'
#' - For a compound like water, it is better to use effective values for Z and A than Braggs rule.
#'     So, for water it is best to apply Z = 10 and A = 18.0158 in a single call to the Bethe formula.
#'     The Bragg rule by taking the Sel/rho(H2O)  = w.H * Sel/rho(H) + w.O * Sel/rho(O) can easily be
#'     a some percent off.
#'
#' - The remaining deviations between the Bethe estimate using Z=2x1+8=10 and A=18.0158 originates from
#'     the delta value used for the density effect.
#'

#' @export
############################################################################
# Water stopping power demonstration
############################################################################
demo.Bragg.rule.test <- function(){
# Created: July 29, 2022
# Revised: July 31, 2022
# Name:    Claus E. Andersen
print("demo.Bragg.rule.test (water at 100 keV, 1 MeV and 10 MeV)")
print("This script should be run manually, line by line.")
print("We compare computations of the mass electronic stopping power")
print("using Bragg's rule (and the MSPs for H and O) and using a")
print("single call to the Bethe formula using effective values for Z = 10 and A = 18.058")
print("We test at 100 keV where there is only a small density effect and at 1 MeV and 10 MeV.")
print("We use the Sternheimer fitted parameters for the density effect in H, O, and H2O.")

dat.H <- list(
  I=19.2,
  Z=1,
  A=1.0079,
  param.Sternheimer.tab.id = "hydrogen-1",
  param.note="Sternheimer et. al 1984, I = 19.2 rho = 8.3748E-5",
  param.C = -9.5835,
  param.X0 = 1.8639,
  param.X1 = 3.2718,
  param.a  = 0.14092,
  param.m  = 5.7273,
  param.delta.X0 = 0.00
)

dat.O <- list(
  I=95,
  Z=8,
  A=16.0,
  param.Sternheimer.tab.id = "oxygen",
  param.note="Sternheimer et. al 1984, I = 95 and rho = 1.3315e-3",
  param.C = -10.7004,
  param.X0 = 1.7541,
  param.X1 = 4.3213,
  param.a  = 0.11778,
  param.m  = 3.2913,
  param.delta.X0 = 0.00
)

dat.H2O <- list(
  I=75,
  Z=10,
  A=18.0158,
  param.Sternheimer.tab.id = "water-1",
  param.note="Sternheimer et. al 1984, I = 75 and rho = 1.000 (liquid)",
  param.C = -3.5017,
  param.X0 = 0.2400,
  param.X1 = 2.8004,
  param.a  = 0.09116,
  param.m  = 3.4773,
  param.delta.X0 =  0.097
)

###############################
# 100 keV
###############################
print('Energy = 100 keV')
MeV <- 0.1

dat.H <- Sternheimer.delta.param(MeV,dat.H)
MSP.H <- electronic.MSP(MeV, dat.H, delta=dat.H$param.delta)


dat.O <- Sternheimer.delta.param(MeV,dat.O)
MSP.O <- electronic.MSP(MeV, dat.O, delta=dat.O$param.delta)

MSP.H2O.Bragg <- MSP.H * 0.11189 + MSP.O * 0.88811
print('Result should be:  4.162491')

print('Set the NIST reference value:')
MSP.H2O.NIST <- 4.115

print('Deviation between the Bragg-rule esimate: 1.15%:')
(MSP.H2O.Bragg - MSP.H2O.NIST)/MSP.H2O.NIST*100

print('We could compute an average density effect and apply it to both H and O.')
print('However this does not make a big difference.')
dat.H2O <- Sternheimer.delta.param(MeV, dat.H2O)
print('Result should be: 0.01380142')

#delta.water <-   0.0
MSP.H2O.Bethe <- electronic.MSP(MeV, dat.H2O, delta=dat.H2O$param.delta)
print('Result should be: 4.11128')

print('Deviation between the Bethe estimate using Z=2x1+8=10 and A=18.0158: -0.09%:')
print('# Note that the delta correction is of little importance for low energies.')
(MSP.H2O.Bethe - MSP.H2O.NIST)/MSP.H2O.NIST*100


###############################
# 1 MeV
###############################
print('Energy = 1 MeV')
MeV <- 1
dat.H <- Sternheimer.delta.param(MeV,dat.H)
MSP.H <- electronic.MSP(MeV, dat.H, delta=dat.H$param.delta)


dat.O <- Sternheimer.delta.param(MeV,dat.O)
MSP.O <- electronic.MSP(MeV, dat.O, delta=dat.O$param.delta)

MSP.H2O.Bragg <- MSP.H * 0.11189 + MSP.O * 0.88811
print('Result should be:  1.8885')

print('Set the NIST reference value:')
MSP.H2O.NIST <- 1.849


print('Deviation between the Bragg-rule esimate: 2.1%:')
(MSP.H2O.Bragg - MSP.H2O.NIST)/MSP.H2O.NIST*100

print('We could compute an average density effect and apply it to both H and O.')
print('However this does not make a big difference.')

dat.H2O <- Sternheimer.delta.param(MeV, dat.H2O)

print('Result = 0.3395958')
print('NIST value = 0.2428')
print('So, no impact of selecting the plasma model (We get the same results as before).')

#delta.water <-   0.0
MSP.H2O.Bethe <- electronic.MSP(MeV, dat.H2O, delta=dat.H2O$param.delta)
MSP.H2O.Bethe
print('Result should be: 1.839749')

print('Deviation between the Bethe esimate using Z=2x1+8=10 and A=18.0158: -0.5%:')
(MSP.H2O.Bethe - MSP.H2O.NIST)/MSP.H2O.NIST*100

print('If we apply the NIST value for delta (0.2428), we get even better agreement:')
MSP.H2O.Bethe.NIST.delta <- electronic.MSP(MeV, dat.H2O, delta=0.2428)

MSP.H2O.Bethe.NIST.delta
print('Result should be: 1.849064')

print('Deviation between the Bethe esimate using Z=2x1+8=10 and A=18.0158: 0.003%:')
(MSP.H2O.Bethe.NIST.delta - MSP.H2O.NIST)/MSP.H2O.NIST*100



###############################
# 10 MeV
###############################
print('Energy = 10 MeV')
MeV <- 10
dat.H <- Sternheimer.delta.param(MeV,dat.H)
MSP.H <- electronic.MSP(MeV, dat.H, delta=dat.H$param.delta)


dat.O <- Sternheimer.delta.param(MeV,dat.O)
MSP.O <- electronic.MSP(MeV, dat.O, delta=dat.O$param.delta)

MSP.H2O.Bragg <- MSP.H * 0.11189 + MSP.O * 0.88811
print('Result should be:  2.237801')

print('Set the NIST reference value:')
MSP.H2O.NIST <- 1.968
# delta.water =  2.906405
# NIST delta =   2.992E+00

print('Deviation between the Bragg-rule esimate: 13.7%:')
(MSP.H2O.Bragg - MSP.H2O.NIST)/MSP.H2O.NIST*100

print('Let us could compute an average density effect and apply it to both H and O.')

dat.H2O <- Sternheimer.delta.param(MeV, dat.H2O)

print('Result =  2.90640483377496')
print('NIST value = 2.992')

#delta.water <-   0.0
MSP.H2O.Bethe <- electronic.MSP(MeV, dat.H2O, delta=dat.H2O$param.delta)
MSP.H2O.Bethe
print('Result: 1.975271')

print('Deviation between the Bethe esimate using Z=2x1+8=10 and A=18.0158: 0.37%:')
(MSP.H2O.Bethe - MSP.H2O.NIST)/MSP.H2O.NIST*100

print('If we apply the NIST value for delta (2.992E+00), we get even better agreement:')
MSP.H2O.Bethe.NIST.delta <- electronic.MSP(MeV, dat.H2O, delta=2.992)

MSP.H2O.Bethe.NIST.delta
print('Result = 1.967959')

print('Deviation between the Bethe esimate using Z=2x1+8=10 and A=18.0158: 0.002%:')
(MSP.H2O.Bethe.NIST.delta - MSP.H2O.NIST)/MSP.H2O.NIST*100

print('Conclusions:')
print('(1) For a compound like water, it is better to use effective values for Z and A than Braggs rule.')
print('    So, for water it is best to apply Z = 10 and A = 18.0158 in a single call to the Bethe formula.')
print('    The Bragg rule by taking the Sel/rho(H2O)  = w.H * Sel/rho(H) + w.O * Sel/rho(O) can easily be')
print('    a few percent off.')
print('(2) The remaining deviations between the Bethe esimate using Z=2x1+8=10 and A=18.0158 originate from')
print('    the delta values used for the density effect.')
} # End demo function for water



