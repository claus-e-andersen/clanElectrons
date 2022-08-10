# Created: July 28, 2022
# Revised: July 31, 2022
# Name   : Claus E. Andersen

# Mass electronic stopping Power for electrons
# Restricted mass electronic stopping Power for electrons
# References:
#  ICRU-90 sec. 4.1.1
#  Andero et al. (2017) p. 74

# The main difference between ICRU-90 values and these functions are in
# the selection of parameters and the computation of the density effect.
# This package uses the simplified approach provided by Sternheimer et al. (1984)
# whereas ICRU-90 uses a more accurate approach. The agreement seems to be
# within 0.5%. The functions in this package can handle the
# density effect corrections in three ways:
#
# (1) We can supply a value for delta directly using the delta.fixed argument.
#     Setting delta.fixed = 1.3 forces the stopping power computations to use this value
#     regardless of all other settings.
# (2) The Sternheimer coefficients are specified for the given material using the
#     Sternheimer et al. (1984) model. If the parameters for the model
#     has already been implemented in the read.Sternheimer.data function, then we get
#     the density effect correction automatically from the the specified Sternheimer.tab.id.
# (3) The Sternheimer coefficients include an 'average' C-value. However, we can also compute this
#     based on the a plasma electron-gas model thereby including an explicit variation with I and rho.
#     If we set C.model to "plasma", we use this option,
#     otherwise we use the 'average' value from Sternheimer directly.
# (4) Finally, we can supply our own Sternheimer coefficients. This is useful if we need
#     computations for materials not already included in the read.Sternheimer.data function.

# Main functions:
#   electronic.MSP.Bethe()             # Compute mass electronic stopping power for electrons
#   restricted.electronic.MSP.Bethe()  # Compute restricted mass electronic stopping power for electrons
#   read.Sternheimer.data()            # Read small database with Sternheimer data
#   delta.Sternheimer()                # Compute delta = the density effect correction
#
# Applications and demonstrations:
#   sensitivity.electronic.MSP.Bethe.dlog.dlogI()
#   demo.stopping.power.validation()
#   demo.stopping.power.plot()
#   demo.stopping.power.for.water()
#   demo.stopping.power.computations()

#' electronic.MSP.Bethe
#' Not so easy as Sternheimer
#' @export
############################################################################
# Mass electronic stopping power
############################################################################
electronic.MSP.Bethe <- function(MeV = 1, dat = NULL(), delta=NA){

#                                 I = 81, Z = 6, A = 12.011, rho = 2.265,
#                                 Sternheimer.tab.id = "graphite-1",
#                                 C.model = "fixed",
#                                 delta.fixed = NA,
#                                 df.Sternheimer = NULL,
#                                 verbose = FALSE){
# Created: July 29, 2022
# Revised: July 31, 2022
# Name:    Claus E. Andersen
# Input:
#   MeV = kinetic energy of electron in MeV (this can be a vector)
#   I   = ionization energy in eV
#   Z   = atomic number (e.g. 6 for graphite)
#   A   = atomic mass (e.g. 12.011 for graphite)
#   rho = density in g/cm3
#   Sternheimer.tab.id = the identified needed to find your material in the database.
#   C.model = "plasma" or "fixed" (the model used by the Sternheimer function).
#   delta.fixed = force this value for the density correction (overwrite Sternheimer).
#   df.Sternheimer = supply Sternheimer data on the fly (if not included in database)
# Output:
#   The mass electronic stopping power for electrons in units of MeV pr. g/cm2
# Sample call:
#   electronic.MSP.Bethe(MeV=1)
#   electronic.MSP.Bethe(MeV=c(0.01,0.1,1,10))
#   electronic.MSP.Bethe(MeV=100,C.model="fixed")
#   electronic.MSP.Bethe(MeV=100,C.model="plasma")
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

# Density effect:
#delta <- delta.Sternheimer(MeV=MeV,I=I, Z=Z, A=A, rho=rho,
#                           Sternheimer.tab.id = Sternheimer.tab.id,
#                           C.model=C.model,
#                           df.Sternheimer=df.Sternheimer,
#                           verbose=verbose)

#if(!is.na(delta.fixed)){
#  # Overwrite Sternheimer
#  if(verbose){
#    print("Message from electronic.MSP.Bethe")
#    print("A delta.fixed value has been supplied, and this will be used in the MSP computations.")
#    print(paste("delta =", delta.fixed))
#  }
#  delta <- delta.fixed
#}

Fminus <- (1-beta^2)*(1 + tau^2/8 - (2*tau+1)*log(2))

P1 <- 2*pi*re^2*E0/beta^2 * Z * Na /  A
P2 <- 2*log(T*1e6/I) + log(1+tau/2) + Fminus - delta

P1 * P2 /( 1E6 * e) * 10000 # MeV pr. g/cm2
} # End function

#' restricted.electronic.MSP.Bethe
#' Not so easy as Sternheimer
#' @export

############################################################################
# Restricted Mass electronic stopping power
############################################################################
restricted.electronic.MSP.Bethe <- function(MeV = 1, delta.keV = 10,
                                            dat=NULL, delta=NA){
#                                            I = 81, Z = 6, A = 12.011, rho = 2.265,
#                                            Sternheimer.tab.id = "graphite-1",
#                                            C.model = "fixed",
#                                            delta.fixed = NA,
#                                            df.Sternheimer = NULL,
#                                           verbose=FALSE){
# Created: July 29, 2022
# Revised: July 31, 2022
# Name:    Claus E. Andersen
# Input:
#   MeV = kinetic energy of electron in MeV (this can be a vector)
#   I   = ionization energy in eV
#   delta.kV = threshold energy in keV for the restriction
#   Z   = atomic number (e.g. 6 for graphite)
#   A   = atomic mass (e.g. 12.011 for graphite)
#   rho = density in g/cm3
#   Sternheimer.tab.id = the identified needed to find your material in the database.
#   C.model = "plasma" or "fixed" (the model used by the Sternheimer function).
#   delta.fixed = force this value for the density correction (overwrite Sternheimer).
#   df.Sternheimer = supply Sternheimer data on the fly (if not included in database)
# Output:
#   The restricted mass electronic stopping power for electrons in units of MeV pr. g/cm2
# Sample call:
#   restricted.electronic.MSP.Bethe(MeV=1,0.1)
#   restricted.electronic.MSP.Bethe(MeV=1,1)
#   restricted.electronic.MSP.Bethe(MeV=1,10)
#   restricted.electronic.MSP.Bethe(MeV=1,100)
#   restricted.electronic.MSP.Bethe(MeV=1,1000)
#   electronic.MSP.Bethe(MeV=1)

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

## Density effect:
#delta <- delta.Sternheimer(MeV=MeV, I=I, Z=Z, A=A, rho=rho,
#                           Sternheimer.tab.id = Sternheimer.tab.id,
#                           C.model = C.model,
#                           df.Sternheimer = df.Sternheimer,
#                           verbose = verbose)
#if(!is.na(delta.fixed)){
#  # Overwrite Sternheimer
#  if(verbose){
#    print("Message from electronic.MSP.Bethe")
#    print("A delta.fixed value has been supplied, and this will be used in the MSP computations.")
#    print(paste("delta =", delta.fixed))
#  }
#  delta <- delta.fixed
#}

# ICRU-90 eq. 4.12
Hminus <-  -1 - beta^2 + log( 4*(1-eta)*eta) + (1-eta)^-1 + (1-beta^2)*(tau^2*eta^2/2 + (2*tau+1) * log(1-eta))

P1 <- 2*pi*re^2*E0/beta^2 * Z * Na /  A
P2 <- 2*log(T*1e6/I) + log(1+tau/2) + Hminus - delta
P1 * P2 /( 1E6 * e) * 10000 # MeV pr. g/cm2
} # End function

#' read.Sternheimer.data
#' Not so easy as Sternheimer
#' @export

############################################################################
# Sternheimer model for density effect: Data for different materials
############################################################################

read.Sternheimer.data <- function(print.wanted=FALSE){
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

#' Sternheimer.delta.param
#'
#' @export
#'
############################################################################
# Sternheimer model for density effect
############################################################################
Sternheimer.delta.param <- function(MeV = 1, dat = NA){
# Created: July 29, 2022
# Revised: August 9, 2022
# Name:    Claus E. Andersen
# Input:
#   MeV = kinetic energy of electron in MeV (this can be a vector)
#   I   = ionization energy in eV
#   Z   = atomic number (e.g. 6 for graphite)
#   A   = atomic mass (e.g. 12.011 for graphite)
#   param.C, paran.X0 etc. : Sternheimer model fit parameters
#   # Output:
#   The Sternheimer density effect correction (delta) needed in the Bethe formula.
# Reference: Andero et al. (2017) p. 74

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
#  print("delta 1")
  delta[ok] <- delta1[ok]
}

ok <- X > X1
if(sum(ok)>0){
#  print("delta 2")
  delta[ok] <- delta2[ok]
}

dat$param.MeV <- MeV
dat$param.delta <- delta

# Return the density effect correction, delta.
dat
} # Sternheimer.delta.param


#' demo.stopping.power.validation
#' Not so easy as Sternheimer
#' @export
#'
############################################################################
# Stopping power validation computations
############################################################################
demo.stopping.power.validation <- function(){
# Created: July 29, 2022
# Revised: July 31, 2022
# Name:    Claus E. Andersen
print("This function should be run manually, line by line.")

dat.C <- list(
    Z    = 6,
    A    = 12.011,
    I    = 81,
    exact.rho =  2.265,
    exact.fvec = c(2/6, 2/6, 2/6),
    exact.Evec = c(288.00, 16.59, 11.26),
    exact.plot = !FALSE,
    param.note="Sternheimer et. al 1984, graphite w. I = 78 and rho = 2.265",
    param.C = -2.8680, param.X0 = -0.0178, param.X1 = 2.3415, param.a  = 0.26142, param.m  = 2.8697,
    param.delta.X0 = 0.12
    )
dat.C <- Sternheimer.set.to.conductor(dat.C)
dat <- dat.C

# Note that the param data are for I = 78  whereas the Bethe eq. is evaluated at I = 81
# Sternheimer parameter model
out1 <- Sternheimer.delta.param(0.01,dat) # Result = 0.005147 ICRU-90 ref = 0.002634
out2 <- Sternheimer.delta.param(0.1,dat)  # Result = 0.055967 ICRU-90 ref = 0.04047
out3 <- Sternheimer.delta.param(1.0,dat)  # Result = 0.820704 ICRU-90 ref = 0.7593
out4 <- Sternheimer.delta.param(10,dat)   # Result = 3.460897 ICRU-90 ref = 3.381

# Stern exact  model
out11 <- Sternheimer.delta.exact(0.01,out1) # Result = NA       ICRU-90 ref = 0.002634
out12 <- Sternheimer.delta.exact(0.1,out2)  # Result = NA       ICRU-90 ref = 0.04047
out13 <- Sternheimer.delta.exact(0.42,out3)  # Result = 0.732938 ICRU-90 ref = 0.7593
out13 <- Sternheimer.delta.exact(1.0,out3)  # Result = 0.732938 ICRU-90 ref = 0.7593
out14 <- Sternheimer.delta.exact(10,out4)   # Result = 3.304996 ICRU-90 ref = 3.381

# Unrestricted mass electronic stopping power, delta = Sternheimer parameter model
electronic.MSP.Bethe(MeV=0.010, out11, out11$param.delta) # Result = 19.983 ICRU-90 = 19.99
electronic.MSP.Bethe(MeV=0.100, out12, out12$param.delta) # Result = 3.6500 ICRU-90 = 3.654
electronic.MSP.Bethe(MeV=1.00,  out13, out13$param.delta) # Result = 1.6007 ICRU-90 = 1.606
electronic.MSP.Bethe(MeV=10.00, out14, out14$param.delta) # Result = 1.7232 ICRU-90 = 1.729

# Unrestricted mass electronic stopping power, delta = Sternheimer exact model
electronic.MSP.Bethe(MeV=0.010, out11, out11$exact.delta) # Result = NA      ICRU-90 = 19.99
electronic.MSP.Bethe(MeV=0.100, out12, out12$exact.delta) # Result = NA      ICRU-90 = 3.654
electronic.MSP.Bethe(MeV=1.00,  out13, out13$exact.delta) # Result = 1.60832 ICRU-90 = 1.606
electronic.MSP.Bethe(MeV=10.00, out14, out14$exact.delta) # Result = 1.7352  ICRU-90 = 1.729

# Supply density correction from ICRU-90 using the ICRU-90 value for delta:
electronic.MSP.Bethe(MeV=0.010, out11, delta=0.002634) # Result = 19.99 ICRU-90 = 19.99
electronic.MSP.Bethe(MeV=0.100, out12, delta=0.04047)  # Result = 3.654 ICRU-90 = 3.654
electronic.MSP.Bethe(MeV=1.00,  out13, delta=0.7593)   # Result = 1.606 ICRU-90 = 1.606
electronic.MSP.Bethe(MeV=10.00, out14, delta=3.381)    # Result = 1.729 ICRU-90 = 1.729
# We now have full agreement with the ICRU-90 values. This implies that the
# source of disagreement identified above arises solely from the density effect
# computations.

} # stopping.power.validation.computations

#' demo.stopping.power.plot
#' Not so easy as Sternheimer
#' @export
#'
############################################################################
# Stopping power demonstration
############################################################################
demo.stopping.power.plot <- function(){
# Created: July 29, 2022
# Revised: July 31, 2022
# Name:    Claus E. Andersen
print("This function should be run manually, line by line.")

xx <- seq(-3,2,length=50)
ee <- 10^xx
yy1 <- restricted.electronic.MSP.Bethe(ee,1)
yy2 <- restricted.electronic.MSP.Bethe(ee,10)
yy3 <- restricted.electronic.MSP.Bethe(ee,100)
yy4 <- restricted.electronic.MSP.Bethe(ee,1000)

df <- rbind(
      data.frame(MeV=ee, Sel.rho = yy1, keV.delta=1),
      data.frame(MeV=ee, Sel.rho = yy2, keV.delta=10),
      data.frame(MeV=ee, Sel.rho = yy3, keV.delta=100),
      data.frame(MeV=ee, Sel.rho = yy4, keV.delta=1000)
)

lattice::xyplot(log10(Sel.rho) ~ log10(MeV),
main="Restricted mass electronic stopping power for graphite",
auto.key=list(title="Delta [keV]",columns=4),
groups=keV.delta,
data=df)
}


#' demo.stopping.power.for.water
#' Not so easy as Sternheimer
#' @export
############################################################################
# Water stopping power demonstration
############################################################################
demo.stopping.power.for.water <- function(){
# Created: July 29, 2022
# Revised: July 31, 2022
# Name:    Claus E. Andersen
print("demo.stopping.power.for.water")
print("This script should be run manually, line by line.")
print("We compare computations of the mass electronic stopping power")
print("using Bragg's rule (and the ST-powers for H and O) and using a")
print("single call to the Bethe formula using effective values for Z = 10 and A = 18.058")


###############################
# 100 keV
###############################
MSP.H <- electronic.MSP.Bethe(MeV=0.1, I=19.2, Z=1, A=1.0079, rho = NA * 8.3748e-5,
                                 Sternheimer.tab.id = "hydrogen-1",
                                 C.model="fixed",
                                 #delta.fixed = delta.water,
                                 #df.Sternheimer=NULL,
                                 verbose = FALSE)


MSP.O <- electronic.MSP.Bethe(MeV=0.1, I=95, Z=8, A=16.000, rho = NA * 1.3315e-3,
                                 Sternheimer.tab.id = "oxygen",
                                 C.model="fixed",
                                 #delta.fixed = delta.water,
                                 #df.Sternheimer=NULL,
                                 verbose = FALSE)

MSP.H20.Bragg <- MSP.H * 0.11189 + MSP.O * 0.88811
# Result:  4.162491

MSP.H20.NIST <- 4.115

# Deviation between the Bragg-rule esimate: 1.15%:
(MSP.H20.Bragg - MSP.H20.NIST)/MSP.H20.NIST*100

# We could computate an average density effect and apply it to both H and O.
# However this does not make a big difference.
delta.water <- delta.Sternheimer(MeV=0.1,  I=75, Z=10.000, A=18.0158, rho = 1,
                                 Sternheimer.tab.id = "water-1",
                                 C.model="fixed")
# Result = 0.01380142

#delta.water <-   0.0
MSP.H20.Bethe <- electronic.MSP.Bethe(MeV=0.1, I=75, Z=10.000, A=18.0158, rho = NA * 1.3315e-3,
                                 Sternheimer.tab.id = "water-1",
                                 C.model="fixed",
                                 #delta.fixed = delta.water,
                                 #df.Sternheimer=NULL,
                                 verbose = FALSE)
MSP.H20.Bethe
# Result: 4.11128

# Deviation between the Bethe esimate using Z=2x1+8=10 and A=18.0158: -0.09%:
# Note that the delta correction is of little importance for low energies.
(MSP.H20.Bethe - MSP.H20.NIST)/MSP.H20.NIST*100



###############################
# 1 MeV
###############################
MSP.H <- electronic.MSP.Bethe(MeV=1, I=19.2, Z=1, A=1.0079, rho = NA * 8.3748e-5,
                                 Sternheimer.tab.id = "hydrogen-1",
                                 C.model="fixed",
                                 #delta.fixed = delta.water,
                                 #df.Sternheimer=NULL,
                                 verbose = FALSE)


MSP.O <- electronic.MSP.Bethe(MeV=1, I=95, Z=8, A=16.000, rho = NA * 1.3315e-3,
                                 Sternheimer.tab.id = "oxygen",
                                 C.model="fixed",
                                 #delta.fixed = delta.water,
                                 #df.Sternheimer=NULL,
                                 verbose = FALSE)

MSP.H20.Bragg <- MSP.H * 0.11189 + MSP.O * 0.88811
# Result:  1.8885

MSP.H20.NIST <- 1.849


# Deviation between the Bragg-rule esimate: 2.1%:
(MSP.H20.Bragg - MSP.H20.NIST)/MSP.H20.NIST*100

# We could computate an average density effect and apply it to both H and O.
# However this does not make a big difference.
delta.water <- delta.Sternheimer(MeV=1, I=75, Z=10, A=18.0158, rho=1,
                                 Sternheimer.tab.id = "water-1",
                                 C.model="fixed",
                                 verbose=TRUE)
# Result = 0.3395958
# NIST value = 0.2428

delta.water <- delta.Sternheimer(MeV=1, I=75, Z=10, A=18.0158, rho=1,
                                 Sternheimer.tab.id = "water-1",
                                 C.model="plasma",
                                 verbose=TRUE)
# Result = 0.3395958
# NIST value = 0.2428
# So, no impact of selecting the plasma model (We get the same results as before).

#delta.water <-   0.0
MSP.H20.Bethe <- electronic.MSP.Bethe(MeV=1, I=75, Z=10.000, A=18.0158, rho = 1,
                                 Sternheimer.tab.id = "water-1",
                                 C.model="fixed",
                                 #delta.fixed = delta.water,
                                 #df.Sternheimer=NULL,
                                 verbose = FALSE)
MSP.H20.Bethe
# Result: 1.839749

# Deviation between the Bethe esimate using Z=2x1+8=10 and A=18.0158: -0.5%:
(MSP.H20.Bethe - MSP.H20.NIST)/MSP.H20.NIST*100

# If we apply the NIST value for delta (0.2428), we get even better agreement:
MSP.H20.Bethe.NIST.delta <- electronic.MSP.Bethe(MeV=1, I=75, Z=10.000, A=18.0158, rho = 1,
                                 Sternheimer.tab.id = "water-1",
                                 C.model="fixed",
                                 delta.fixed = 0.2428,
                                 #df.Sternheimer=NULL,
                                 verbose = FALSE)
MSP.H20.Bethe.NIST.delta
# 1.849064

# Deviation between the Bethe esimate using Z=2x1+8=10 and A=18.0158: 0.003%:
(MSP.H20.Bethe.NIST.delta - MSP.H20.NIST)/MSP.H20.NIST*100



###############################
# 10 MeV
###############################
MSP.H <- electronic.MSP.Bethe(MeV=10, I=19.2, Z=1, A=1.0079, rho = NA * 8.3748e-5,
                                 Sternheimer.tab.id = "hydrogen-1",
                                 C.model="fixed",
                                 #delta.fixed = delta.water,
                                 #df.Sternheimer=NULL,
                                 verbose = FALSE)


MSP.O <- electronic.MSP.Bethe(MeV=10, I=95, Z=8, A=16.000, rho = NA * 1.3315e-3,
                                 Sternheimer.tab.id = "oxygen",
                                 C.model="fixed",
                                 #delta.fixed = delta.water,
                                 #df.Sternheimer=NULL,
                                 verbose = FALSE)

MSP.H20.Bragg <- MSP.H * 0.11189 + MSP.O * 0.88811
# Result:  2.237801

MSP.H20.NIST <- 1.968
# delta.water =  2.906405
# NIST delta =   2.992E+00

# Deviation between the Bragg-rule esimate: 13%:
(MSP.H20.Bragg - MSP.H20.NIST)/MSP.H20.NIST*100

# Let's could computate an average density effect and apply it to both H and O.

delta.water <- delta.Sternheimer(MeV=10, I=75, Z=10, A=18.0158, rho=1,
                                 Sternheimer.tab.id = "water-1",
                                 C.model="fixed",
                                 verbose=TRUE)

# Result =  2.90640483377496
# NIST value = 2.992

MSP.H.H20.delta  <- electronic.MSP.Bethe(MeV=10, I=19.2, Z=1, A=1.0079, rho = NA * 8.3748e-5,
                                 Sternheimer.tab.id = "hydrogen-1",
                                 C.model="fixed",
                                 delta.fixed = delta.water,
                                 #df.Sternheimer=NULL,
                                 verbose = FALSE)


MSP.O.H20.delta  <- electronic.MSP.Bethe(MeV=10, I=95, Z=8, A=16.000, rho = NA * 1.3315e-3,
                                 Sternheimer.tab.id = "oxygen",
                                 C.model="fixed",
                                 delta.fixed = delta.water,
                                 #df.Sternheimer=NULL,
                                 verbose = FALSE)

MSP.H20.Bragg.H20.delta <- MSP.H.H20.delta  * 0.11189 + MSP.O.H20.delta  * 0.88811
# Result = 1.98952


#delta.water <-   0.0
MSP.H20.Bethe <- electronic.MSP.Bethe(MeV=10, I=75, Z=10.000, A=18.0158, rho = 1,
                                 Sternheimer.tab.id = "water-1",
                                 C.model="fixed",
                                 #delta.fixed = delta.water,
                                 #df.Sternheimer=NULL,
                                 verbose = FALSE)
MSP.H20.Bethe
# Result: 1.975271

# Deviation between the Bethe esimate using Z=2x1+8=10 and A=18.0158: 0.37%:
(MSP.H20.Bethe - MSP.H20.NIST)/MSP.H20.NIST*100

# If we apply the NIST value for delta (2.992E+00), we get even better agreement:
MSP.H20.Bethe.NIST.delta <- electronic.MSP.Bethe(MeV=10, I=75, Z=10.000, A=18.0158, rho = 1,
                                 Sternheimer.tab.id = "water-1",
                                 C.model="fixed",
                                 delta.fixed = 2.992E+00,
                                 #df.Sternheimer=NULL,
                                 verbose = FALSE)
MSP.H20.Bethe.NIST.delta
# 1.967959

# Deviation between the Bethe esimate using Z=2x1+8=10 and A=18.0158: 0.002%:
(MSP.H20.Bethe.NIST.delta - MSP.H20.NIST)/MSP.H20.NIST*100

# Conclusions
# (1) For a compound like water, it is better to use effective values for Z and A than Bragg'r rult.
#     So, for water it is best to apply Z = 10 and A = 18.0158 in a single call to the Bethe formula.
#     The Bragg rule by taking the Sel/rho(H2O)  = w.H * Sel/rho(H) + w.O * Sel/rho(O) can easily be
#     a few percent off.
# (2) The remaining deviations between the Bethe esimate using Z=2x1+8=10 and A=18.0158 originate from
#     the delta values used for the density effect.

} # End demo function for water




#' demo.stopping.power.computations
#' Not so easy as Sternheimer
#' @export
############################################################################
# Demonstration of misc stopping power computations
############################################################################

demo.stopping.power.computations <- function(){
# Created: July 29, 2022
# Revised: July 31, 2022
# Name:    Claus E. Andersen
print("This function should be run manually, line by line.")

print("demo.stopping.power.computations")
print("This script should be run manually, line by line.")

read.Sternheimer.data()
delta.Sternheimer()
electronic.MSP.Bethe()
electronic.MSP.Bethe(1,C.model="fixed",verbose=TRUE)
electronic.MSP.Bethe(1,C.model="plasma",verbose=TRUE)

# Note that C.model works regardless if the Sternheimer parameters are supplied on the fly or extracted from
# the standard table.

df.Sternheimer.junk <-
data.frame(Z=6, Sternheimer.tab.id.sel ="Junk material", note="", C = -2000.9925, X0 = -0.10351, X1 = 32.486, a  = 0.20240, m  = 3.0036, delta.X0 = 0.10)
electronic.MSP.Bethe(1,C.model="fixed",df.Sternheimer = df.Sternheimer.junk, verbose=TRUE)

df.Sternheimer.junk <-
data.frame(Z=6, Sternheimer.tab.id.sel ="Junk material", note="", C = -2000.9925, X0 = -0.10351, X1 = 32.486, a  = 0.20240, m  = 3.0036, delta.X0 = 0.10)
electronic.MSP.Bethe(1,C.model="plasma",df.Sternheimer = df.Sternheimer.junk, verbose=TRUE)

# Note that we will always use delta.fixed value (if supplied) - regardless of all other parameters.
electronic.MSP.Bethe(1,C.model="fixed", df.Sternheimer = df.Sternheimer.junk, verbose=TRUE, delta.fixed=0)
electronic.MSP.Bethe(1,C.model="plasma",df.Sternheimer = df.Sternheimer.junk, verbose=TRUE, delta.fixed=0)



restricted.electronic.MSP.Bethe()

## delta.Sternheimer(Sternheimer.tab.id = "Unknown material")

restricted.electronic.MSP.Bethe(c(0.01,0.1,1,10),10000, df.Sternheimer=df.Sternheimer.junk)


electronic.MSP.Bethe(c(0.1,1,10))
electronic.MSP.Bethe(c(0.1,1,10),df.Sternheimer=df.Sternheimer.junk)
}#

