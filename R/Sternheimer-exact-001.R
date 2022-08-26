#' @title clanElectrons
#' @description  Package for computation of electronic stopping power and density effect (Sternheimer)
#' Dummy function to ease acces to the help index.
#'
#'  clanElectrons()
#'
#'  ?clanElectrons # gives you an index of functions in package
#'
#'  #' Status (August 10, 2022):
#'
#'  1. Bethe formula and Sternheimer.exact.delta agres perfectly with ICRU-90
#'  for water, but NOT for graphite.
#'
#'  2. More work needed for the lower energies for insulators and conductors.
#'  Clean up of examples are needed.
#'
#'  3. Consider if material data could be stored.
#'
#'  4. Perhaps we should make a set of functions:
#'
#'   demo.MSP.water
#'
#'   demo.MSP.graphite
#'
#'   demo.alanine
#'
#'   demo.MSP.Al
#'
#'   demo.MSP.Be
#'
#'   etc...ß
#'
#'
#' @export
clanElectrons <- function(){
# Dummy function
}

#' @title Sternheimer.set.to.conductor
#' @description  Ascertain that material is treated as conductor in Sternheimer delta computations
#'
#'  dat <- Sternheimer.set.to.conductor(dat)
#'
#' The binding energy of the outer subshell is set to zero
#' @export

#' @title Sternheimer.f.root.mu.st <-
#' @description  Helper function (mu.st) for exact Sternheimer density correction
#'
#' # Equation suitable for finding mu.st (eq. 4.29 in ICRU-90)
#' @export
Sternheimer.f.root.mu.st <- function(mu.st,dat){
# Helper function (mu.st) for exact Sternheimer density correction
# Created: Aug 7, 2022
# Revised: Aug 25, 2022
# Name   : Claus E. Andersen
# Modified version of eq. 4.29 in ICRU-90
# suitable for root finding: The requested
# x = mu.st will fullfill: f.root.mu.st(x) = 0
fvec <- dat$fvec
Evec <- dat$Evec
Ep   <- dat$Ep
nlev <- dat$nlev
nc   <- dat$nc
fnc  <- nc/dat$Z
I    <- dat$I
ans <- 0
fn <- fvec[nlev]
for(i in 1:(nlev)){
  if(fvec[i]>0){
    ans <- ans + fvec[i] * log(    ( (mu.st * Evec[i])^2 + 2/3 * fvec[i] * Ep^2 )^0.5   )
  }# if
}# loop
#ans <- ans + fn * log((fn)^0.5 * Ep)
if(nc >0){ ans <- ans + fnc * log(fnc^0.5*Ep)  }
ans <- ans - log(I)
ans
}# Sternheimer.f.root.mu.st

###########################################################################

#' @title Sternheimer.f.root.L
#' @description Helper function (L) for exact Sternheimer density correction
#' # Equation suitable for finding L (ell) (eq. 4.28 in ICRU-90)
#' @export
Sternheimer.f.root.L <- function(L,dat){
# Helper function (L) for exact Sternheimer density correction
# Created: Aug 7, 2022
# Revised: Aug 25, 2022
# Name   : Claus E. Andersen
# Modified version of eq. 4.28 in ICRU-90
# suitable for root finding: The requested
# x = L will fullfill: f.root.L(x) = 0
fvec  <- dat$fvec
Evec  <- dat$Evec
mu.st <- dat$mu.st
Ep    <- dat$Ep
beta  <- dat$beta
nlev  <- dat$nlev
nc    <- dat$nc
fnc   <- nc/dat$Z

ans <- 0
for(i in 1:(nlev)){
  if(fvec[i]>0){
  ans <- ans + fvec[i] / ( (mu.st*Evec[i]/Ep)^2 + L^2 )
  }# if
}# loop
if(nc>0){ ans <- ans + fnc / L^2 }
ans <- ans + 1 - 1/beta^2
ans
}# Sternheimer.f.root.L


###############################################################
###############################################################


#' @title Sternheimer.delta.exact
#' @description  Computation of exact Sternheimer density correction.
#'
#' By exact, we just mean that this is not the simple
#' parameter fitting solution. The ICRU-90 terminology
#' (p. 26) is followed as close as possible).
#'
#' @param dat list with parameters
#' @param dat$MeV = electron kinetic energy
#' @param dat$nlev = number of sub-shells
#' @param dat$Z = atomic number
#' @param dat$A = atomic mass
#' @param dat$rho.density = density in g/cm3
#' @param dat$fvec = sub-shell occupancy level (used in computation)
#' @param dat$Evec = sub-shell binding energy (used in computation)
#' @param dat$fvec.org = sub-shell occupancy level (NOT used in computation)
#' @param dat$Evec.org = sub-shell binding energy (NOT used in computation)
#' @param dat$I = mean excitation energy in eV
#' @param dat$plot.wanted = TRUE or FALSE,
#' @details Notes:
#' (1) The equations for mu.st and L are solved using a search strategy
#'     (coarse + fine + regression) which may fail, for example, if the first interval
#'     dose not include the root. Graphs can ne plotted to inspect  the procedure.
#' (2) MeV is the kinetic energy and only one value can be computeted (not vectorized).
#' (3) Only computations for single atoms has been implemented (not multi-atom materials).
#' (4) The data fit well the ICRU-37 delta values for high energies (e.g. at 100 MeV), but
#'     less well for lower energies (e.g. at 1 MeV). The source of deviation has not been
#'     identified.
#' (5) The main use of this work is to try to gain some insights into the procedure.
#'     Like how conductors are handled. And that insulators have delta = 0 for energies
#'     below a certain threshold beta0 (not roots).
#' (6) Why is this insight relevant? This is needed as the text is not super clear and
#'     since the  whole concept of dialetric / optical modelling is new to me (e.g. the
#'     use of plasma energy).
#'
#' Main references:
#' (1) ICRU-90: Key data for ionizing-radiation dosimetry: Measurement standards
#'     and applications (2014/2016).
#' (2) DENSITY EFFECT FOR THE IONIZATION LOSS OF CHARGED PARTICLES
#'     IN VARIOUS SUBSTANCES
#'     by R. M. STERNHEIMER (Brookhaven), M. J. BERGER (NBS) and S. M. SELTZER (NBS).
#'     ATOMIC DATA AND NUCLEAR DATA TABLES 30,26 l-27 1 ( 1984)
#' (3) G4DensityEffectCalculator.cc by Matthew Strait <straitm@umn.edu> 2019
#' @export
Sternheimer.delta.exact <- function(MeV=1, dat=NULL, mu.solver.parm = NULL,  L.solver.parm = NULL) {
# Created: Aug 7, 2022
# Revised: Aug 26, 2022
# Name   : Claus E. Andersen
#
# Computation of exact Sternheimer density correction (delta).
#
# By exact, we just mean that this is not the simple
# parameter fitting solution. The ICRU-90 terminology
# (p. 26) is followed as close as possible).

# Notes:
# (1) The equations for mu.st and L are solved using a search strategy
#     (coarse + fine + regression) which may fail, for example, if the first interval
#     dose not include the root. Graphs can ne plotted to inspect  the procedure.
# (2) MeV is the kinetic energy and only one value can be computeted (not vectorized).
# (3) Only computations for single atoms has been implemented (not multi-atom materials).
# (4) The data fit well the ICRU-37 delta values for high energies (e.g. at 100 MeV), but
#     less well for lower energies (e.g. at 1 MeV). The source of deviation has not been
#     identified.
# (5) The main use of this work is to try to gain some insights into the procedure.
#     Like how conductors are handled. And that insulators have delta = 0 for energies
#     below a certain threshold beta0 (not roots).
# (6) Why is this insight relevant? This is needed as the text is not super clear and
#     since the  whole concept of dialetric / optical modelling is new to me (e.g. the
#     use of plasma energy).
#
# Main references:
# (1) ICRU-90: Key data for ionizing-radiation dosimetry: Measurement standards
#     and applications (2014/2016).
# (2) DENSITY EFFECT FOR THE IONIZATION LOSS OF CHARGED PARTICLES
#     IN VARIOUS SUBSTANCES
#     by R. M. STERNHEIMER (Brookhaven), M. J. BERGER (NBS) and S. M. SELTZER (NBS).
#     ATOMIC DATA AND NUCLEAR DATA TABLES 30,26 l-27 1 ( 1984)
# (3) G4DensityEffectCalculator.cc by Matthew Strait <straitm@umn.edu> 2019


  if(is.null(mu.solver.parm)){
    mu.solver.parm <- list(
      mu.st.min = 0,
      mu.st.max = 20,
      mu.st.N   = 10000,
      mu.st.eps = 6e-05
    )
  }

  if(is.null(L.solver.parm)){
    L.solver.parm <- list(
      L.min = 0.02,
      L.max = 4000,
      L.N   = 80000,
      L.eps  = 1e-03
    )
  }


E0       <- 0.51099895000
dat$beta <- (1 - (E0/(E0+MeV))^2)^0.5
dat$Ep   <- 28.8159 * (dat$exact.rho *dat$Z/dat$A)^0.5
dat$nlev <- length(dat$fvec)
nc       <- dat$nc
fnc      <- nc/dat$Z

dat$exact.MeV <- MeV
dat$mu.st.solution <- "None"


####################################################
# Part 1: Find mu.st (solve eq. 4.29 in ICRU-90)
####################################################
# Search for mu.st (coarse)

xx <- seq(mu.solver.parm$mu.st.min, mu.solver.parm$mu.st.max, length = mu.solver.parm$mu.st.N)
yy <- Sternheimer.f.root.mu.st(xx,dat)

ok.pos <- yy>0
ok.neg <- yy<0

if(sum(ok.pos)> 0 & sum(ok.neg)>0){
# Search for mu.st (fine)
mu.st.min <- max(xx[ok.neg])
mu.st.max <- min(xx[ok.pos]) + 0.001
mu.st.N   <- mu.solver.parm$mu.st.N
mu.st.eps <- mu.solver.parm$mu.st.eps

xx <- seq(mu.st.min, mu.st.max, length = mu.st.N)
yy <- Sternheimer.f.root.mu.st(xx,dat)

ok <- abs(yy) < mu.st.eps & xx>0
mu.st.root <- NA
if(sum(ok)>0){
  # Use regression to find final value (root)
  fm <- lm(xx[ok]~yy[ok])
  mu.st.root <- as.numeric(coefficients(fm)[1])
  dat$mu.st.solution <- "Regression"
}

df <- data.frame(mu.st=xx,val=yy)

plt.mu.st <- lattice::xyplot(val ~ mu.st,
  data=df,
  panel=function(x,y,...){
    lattice::panel.xyplot(x,y,...)
    lattice::panel.abline(h=0,lty="dashed")
    lattice::panel.abline(v=mu.st.root,lty="dashed")
  }
)
if(dat$exact.plot)(print(plt.mu.st))
}# both pos and neg

dat$mu.st <- mu.st.root

####################################################
# PART 2: Find L (solve eq. 4.28 in ICRU-90)
####################################################
# Search interval for L (coarse)
dat$L.solution <- "None"
L.root <- NA

xx <- seq(L.solver.parm$L.min, L.solver.parm$L.max, length = L.solver.parm$L.N)
yy <- Sternheimer.f.root.L(xx,dat)

ok.pos <- yy>0
ok.neg <- yy<0

if(sum(ok.pos)> 0 & sum(ok.neg)>0){
# Search interval for L (fine)
L.min <- max(xx[ok.pos])
L.max <- min(xx[ok.neg])
L.eps <- L.solver.parm$L.eps

xx <- seq(L.min ,L.max, length=L.solver.parm$L.N)
yy <- Sternheimer.f.root.L(xx,dat)

##ok <- abs(yy) < L.eps & xx>0
ok <- abs(yy) < L.eps & yy < 0.5 & yy > -0.5
L.root <- NA
if(sum(ok)>0){
  # Use regression to find final value (root)
  fm <- lm(xx[ok]~yy[ok])
  L.root <- as.numeric(coefficients(fm)[1])
  dat$L.solution <- "Regression"
}


} # both pos and neg

df <- data.frame(L=xx,val=yy)

dat$L  <- L.root

L     <- dat$L
mu.st <- dat$mu.st
Evec  <- dat$Evec
fvec  <- dat$fvec
beta  <- dat$beta
nlev  <- dat$nlev
Ep    <- dat$Ep


plt.L <- lattice::xyplot(val ~ L,
                         data=df,
                         main="Sternheimer.delta.exact, root finding (L equation)",
                         panel=function(x,y,...){
                           lattice::panel.xyplot(x,y,...)
                           lattice::panel.abline(h=0,lty="dashed")
                           lattice::panel.abline(v=L.root,lty="dashed")
                         }
)

if(dat$exact.plot){print(plt.L)}


####################################################
# Part 3: Find delta using eq. 4.27 in ICRU-90
####################################################
Lvec <- (  (mu.st*Evec/Ep)^2 + 2/3*fvec ) ^0.5
if(Evec[nlev]==0){
  Lvec[nlev] <- (fvec[nlev])^0.5
}

dat$Lvec <- Lvec

ans <- 0
for(i in 1:nlev){
  if(fvec[i]>0){
  ans <- ans + fvec[i] * log(1 + L^2/Lvec[i]^2)
}# if
}# loop
# Important (?): The final part in 4.27 is not part of the summation!!!
if(dat$nc>0){ans <- ans + fnc * log(1 + L^2/fnc) }
delta <- ans - L^2 * (1-beta^2)

# Done
dat$exact.delta <- delta

if(dat$exact.plot){clanLattice::print.trellis.plots(list(plt.mu.st,plt.L),1)}

dat
} #Sternheimer.delta.exact.function


#############################################################################


########################################################################################

#' @title demo.Sternheimer.delta.exact.plot
#' Delta vs. energy for two different models (Al: En=0 or not)
#' @export
demo.Sternheimer.delta.exact.plot <- function(){
# Created: July 29, 2022
# Revised: July 31, 2022
# Revised: August 25, 2022
# Name:    Claus E. Andersen
print("This function should be run manually, line by line.")


dat.Al.model0 <- list(
exact.plot = FALSE,
Z    = 13,
A    = 26.98154,
exact.rho =  2.699,
nc = 0,
fvec = c(2/13, 2/13 ,2/13, 2/13, 2/13, 3/13),
Evec = c( 1564.0 , 121.0, 77.0, 77.0, 10.62, 5.986),
I = 166.0)

dat.Al.model1 <- list(
  exact.plot = FALSE,
  Z    = 13,
  A    = 26.98154,
  exact.rho =  2.699,
  nc = 1,
  fvec = c(2/13, 2/13 ,2/13, 2/13, 2/13, 2/13),
  Evec = c( 1564.0 , 121.0, 77.0, 77.0, 10.62, 5.986),
  I = 166.0)

dat.Al.model2 <- list(
  exact.plot = FALSE,
  Z    = 13,
  A    = 26.98154,
  exact.rho =  2.699,
  nc = 2,
  fvec = c(2/13, 2/13 ,2/13, 2/13, 2/13, 1/13),
  Evec = c( 1564.0 , 121.0, 77.0, 77.0, 10.62, 5.986),
  I = 166.0)


dat.Al.model3 <- list(
  exact.plot = FALSE,
  Z    = 13,
  A    = 26.98154,
  exact.rho =  2.699,
  nc = 3,
  fvec = c(2/13, 2/13 ,2/13, 2/13, 2/13, 0/13),
  Evec = c( 1564.0 , 121.0, 77.0, 77.0, 10.62, 5.986),
  I = 166.0)


xx <- seq(-2,3,length=50)
ee <- 10^xx
yy0 <- 0 * xx
yy1 <- 0 * xx
yy2 <- 0 * xx
yy3 <- 0 * xx

dat <- dat.Al.model0
for(ii in 1:length(ee)){
  MeV <- ee[ii]
  dat.out <- Sternheimer.delta.exact(MeV,dat)
  yy0[ii] <- dat.out$exact.delta
}

dat <- dat.Al.model1
for(ii in 1:length(ee)){
  MeV <- ee[ii]
  dat.out <- Sternheimer.delta.exact(MeV,dat)
  yy1[ii] <- dat.out$exact.delta
}

dat <- dat.Al.model2
for(ii in 1:length(ee)){
  MeV <- ee[ii]
  dat.out <- Sternheimer.delta.exact(MeV,dat=dat)
  yy2[ii] <- dat.out$exact.delta
}

dat <- dat.Al.model3
for(ii in 1:length(ee)){
  MeV <- ee[ii]
  dat.out <- Sternheimer.delta.exact(MeV,dat=dat)
  yy3[ii] <- dat.out$exact.delta
}


df0 <- data.frame(MeV=ee, material = "Al (insulator model)", delta=yy0)
df1 <- data.frame(MeV=ee, material = "Al (conductor model, nc = 1)", delta=yy1)
df2 <- data.frame(MeV=ee, material = "Al (conductor model, nc = 2)", delta=yy2)
df3 <- data.frame(MeV=ee, material = "Al (conductor model, nc = 3)", delta=yy3)
df <- rbind(df0,df1,df2,df3)

plt <- lattice::xyplot(delta ~ log10(MeV) | material,
par.strip.text=list(cex=1.5),
main="Density correction (delta) from exact Sternheimer procedure\n Red points are from ICRU-37. Aluminium (Z=13).",
panel=function(x,y,...){
lattice::panel.xyplot(x,y,...)
xx <- c(0.1,0.4,1,4,10,40,100,400,1000)
yy <- c(0.01513, 0.1190,0.3339,1.183,2.384,4.669,6.363,9.091,10.92)
lattice::panel.points(log10(xx),yy,col="red",pch=16,cex=0.6)
},
data=df)

print(plt)

plt2 <- lattice::xyplot(delta ~ log10(MeV), groups= material,
                        auto.key=list(columns=2),
                       par.strip.text=list(cex=1.5),
                       main="Density correction (delta) from exact Sternheimer procedure\n Red points are from ICRU-37. Aluminium (Z=13).",
                       panel=function(x,y,...){
                         lattice::panel.xyplot(x,y,...)
                         xx <- c(0.1,0.4,1,4,10,40,100,400,1000)
                         yy <- c(0.01513, 0.1190,0.3339,1.183,2.384,4.669,6.363,9.091,10.92)
                         lattice::panel.points(log10(xx),yy,col="red",pch=16,cex=0.6)
                       },
                       data=df)

print(plt2)
} # demo.Sternheimer.delta.exact.plot



#' @title demo.Sternheimer.water
#' Computation of electronic stopping power and density effect for liqud water
#' using Sternheimer model as described in ICRU-90
#' @export
demo.Sternheimer.water <- function(){
# Created: August 9, 2022
# Revised: August 25, 2022
# Name:    Claus E. Andersen

# How to compute the density-correction for a compound
# like water? Add the Z for all atoms involved.

# Arrange the fvec and the Evec atom by atom.
# Compute electron subshell occopancy factor as the number of electrons
# divided by the total Z. So, for water we have Z.sum = 10 electrons.
# First we consider the two hydrogen atoms, they have 13.6 eV binding energy
# and we therefore set fvec[1] to 2/10 and Evec[1] to 13.6. Then we have the
# 8 electrons in oxygen: fvec[2] = 2, fvec[3] = 2, and fvec[4] = 4, with
# binding energies:  Evec[2] = 538.0, Evec[3] = 28.48, and Evec[4] = 13.62.
# Compounds should be treated as an insulator.

# In this example, we use the recommended values of I = 78 eV for water
# and we also set the density to 0.998 g/cm3 which has some implications
# for the density correction correction.

dat.H2O <- list(
  Z    = 10,
  A    = 18.0158,
  I    = 78,
  exact.rho =  0.998,
  fvec = c(2/10, 2/10, 2/10, 4/10),
  Evec = c(13.6, 538.0, 28.48, 13.62),
  nc = 0,  # Always set compounds to insulators
  exact.plot = FALSE,
  param.note="Sternheimer et. al 1984, water (liquid) I = 75 and rho = 1.000 ",
  param.C = -3.5017, param.X0 = 0.2400, param.X1 = 2.8004, param.a  = 0.09116, param.m  = 3.4773,
  param.delta.X0 = 0.097
)

dat <- dat.H2O
############################
# 800 keV
############################
MeV <- 0.8
dat <- Sternheimer.delta.exact(MeV, dat)
xx0 <- electronic.MSP.Bethe(MeV, dat, delta = 0) # No density effect correction
xx <- electronic.MSP.Bethe(MeV, dat, delta = dat$exact.delta) # Exact Sternheimer density correction
df1 <- data.frame(MeV = MeV, I.eV = dat$I, rho=dat$exact.rho,  MSP.R0 = xx0,
                  MSP.R = xx, MSP.ICRU90=1.880,
                  delta.R=dat$exact.delta, delta.ICRU90=0.1005,
                  mu.st = dat$mu.st,
                  L = dat$L)

############################
# 1 MeV
############################
MeV <- 1.0
dat <- Sternheimer.delta.exact(MeV, dat)
xx0 <- electronic.MSP.Bethe(MeV, dat, delta = 0) # No density effect correction
xx  <- electronic.MSP.Bethe(MeV, dat, delta = dat$exact.delta) # Exact Sternheimer density correction
df2 <- data.frame(MeV = MeV, I.eV = dat$I, rho=dat$exact.rho,  MSP.R0 = xx0,
                  MSP.R = xx, MSP.ICRU90 = 1.845,
                  delta.R = dat$exact.delta, delta.ICRU90 = 0.2086,
                  mu.st = dat$mu.st,
                  L = dat$L)

############################
# 10 MeV
############################
MeV <- 10.0
dat <- Sternheimer.delta.exact(MeV, dat)
xx0 <- electronic.MSP.Bethe(MeV, dat, delta = 0) # No density effect correction
xx  <- electronic.MSP.Bethe(MeV, dat, delta = dat$exact.delta) # Exact Sternheimer density correction
df3 <- data.frame(MeV = MeV, I.eV = dat$I, rho=dat$exact.rho,  MSP.R0 = xx0,
                  MSP.R = xx, MSP.ICRU90 = 1.967,
                  delta.R = dat$exact.delta, delta.ICRU90 = 2.928,
                  mu.st = dat$mu.st,
                  L = dat$L)

############################
# 100 MeV
############################
MeV <- 100.0
dat <- Sternheimer.delta.exact(MeV, dat)
xx0 <- electronic.MSP.Bethe(MeV, dat, delta = 0) # No density effect correction
xx  <- electronic.MSP.Bethe(MeV, dat, delta = dat$exact.delta) # Exact Sternheimer density correction
df4 <- data.frame(MeV = MeV, I.eV = dat$I, rho=dat$exact.rho,  MSP.R0 = xx0,
                  MSP.R = xx, MSP.ICRU90 = 2.202,
                  delta.R = dat$exact.delta, delta.ICRU90 = 6.998,
                  mu.st = dat$mu.st,
                  L = dat$L)


############################
# 1000 MeV
############################
MeV <- 1000.0
dat <- Sternheimer.delta.exact(MeV, dat)
xx0 <- electronic.MSP.Bethe(MeV, dat, delta = 0) # No density effect correction
xx  <- electronic.MSP.Bethe(MeV, dat, delta = dat$exact.delta) # Exact Sternheimer density correction
df5 <- data.frame(MeV = MeV, I.eV = dat$I, rho=dat$exact.rho,  MSP.R0 = xx0,
                  MSP.R = xx, MSP.ICRU90 = 2.401,
                  delta.R = dat$exact.delta, delta.ICRU90 = 11.58,
                  mu.st = dat$mu.st,
                  L = dat$L)

# Combine all results:
df <- rbind(df1,df2,df3,df4,df5)

############################################################################
# Main results
############################################################################

# Excellent agreement with between the clanElectrons computations
# and ICRU-90 values for water. MSP.R0 is the mass electronic
# stopping power without density effect correction (delta=0).
#
#     MeV I.eV   rho   MSP.R0    MSP.R MSP.ICRU90    delta.R delta.ICRU90
# 1 8e-01   78 0.998 1.890531 1.880437      1.880  0.1004488       0.1005
# 2 1e+00   78 0.998 1.864880 1.844806      1.845  0.2086075       0.2086
# 3 1e+01   78 0.998 2.216852 1.966726      1.967  2.9279906       2.9280
# 4 1e+02   78 0.998 2.798672 2.202281      2.202  6.9977588       6.9980
# 5 1e+03   78 0.998 3.387159 2.400502      2.401 11.5772526      11.5800


############################################################################
# Alternative computations.
############################################################################

# Computation with density set to 1 (ICRU-90 values are for 78/0.998)
#    MeV I.eV rho   MSP.R0    MSP.R MSP.ICRU90    delta.R delta.ICRU90
#1 8e-01   78   1 1.890531 1.880374      1.880  0.1010731       0.1005
#2 1e+00   78   1 1.864880 1.844726      1.845  0.2094346       0.2086
#3 1e+01   78   1 2.216852 1.966590      1.967  2.9295878       2.9280
#4 1e+02   78   1 2.798672 2.202113      2.202  6.9997265       6.9980
#5 1e+03   78   1 3.387159 2.400331      2.401 11.5792542       11.580

# Computation with I et to 75 (ICRU-90 values are for 78/0.998)
#    MeV I.eV   rho   MSP.R0    MSP.R MSP.ICRU90    delta.R delta.ICRU90
#1 8e-01   75 0.998 1.898414 1.885749      1.880  0.1260246       0.1005
#2 1e+00   75 0.998 1.872428 1.849151      1.845  0.2418950       0.2086
#3 1e+01   75 0.998 2.223553 1.968077      1.967  2.9906215       2.9280
#4 1e+02   75 0.998 2.805357 2.202391      2.202  7.0749039       6.9980
#5 1e+03   75 0.998 3.393844 2.400503      2.401 11.6556801       11.580

# Computation with I et to 75 and rho =1 (ICRU-90 values are for 78/0.998)
# These results should align with ICRU-37.
#    MeV I.eV rho   MSP.R0    MSP.R MSP.ICRU90    delta.R delta.ICRU90
#1 8e-01   75   1 1.898414 1.885681      1.880  0.1267066       0.1005
#2 1e+00   75   1 1.872428 1.849067      1.845  0.2427676       0.2086
#3 1e+01   75   1 2.223553 1.967940      1.967  2.9922215       2.9280
#4 1e+02   75   1 2.805357 2.202223      2.202  7.0768740       6.9980
#5 1e+03   75   1 3.393844 2.400333      2.401 11.6576817       11.580

df
} # water (demo)


#' @title demo.Sternheimer.graphite
#' Computation of electronic stopping power and density effect for graphite
#' using Sternheimer model as described in ICRU-90
#' @export
demo.Sternheimer.graphite <- function(){
  # Created: August 25, 2022
  # Revised: August 25, 2022
  # Name:    Claus E. Andersen

  # How to compute the density-correction for a conductor
  # like graphite? Move one or more electrons to nc (i.e.
  # binding energy zero).

  # Compounds should be treated as an insulator.

  # In this example, we use the recommended values of I = 81 eV for graphite
  # and we also set the density to 2.265 g/cm3 which has some implications
  # for the density correction correction.

  dat.graphite <- list(
    Z    = 6,       # Atomic number
    A    = 12.011,  # Atomic mass
    I    = 81,      #78,       # Mean excitation energy in eV
    nc = 1,
    exact.rho =  2.265,         # Density in g/cm3 only needed for the exact density-effect correction.
    fvec = c(2/6, 2/6,1/6),     # Occupation fractions for the subshells in C
    Evec = c(288, 16.59,11.26), # Binding energies of subshells from Carlson (1975), see ICRU-90.
    exact.plot = FALSE          # Supplementary plots related to the root finding in the exact density correction
  )

  dat <- dat.graphite
  ############################
  # 1 keV
  ############################
  MeV <- 0.001 # Electron kinetic energy
  dat <- Sternheimer.delta.exact(MeV, dat)
  xx0 <- electronic.MSP.Bethe(MeV, dat, delta = 0)               # Compute MSP without density effect correction (delta = 0)
  xx  <- electronic.MSP.Bethe(MeV, dat, delta = dat$exact.delta) # Compute MSP with exact Sternheimer density correction

  df100 <- data.frame(
    MeV = MeV,
    I.eV = dat$I,
    rho=dat$exact.rho,
    MSP.R0 = xx0,
    MSP.R = xx,
    MSP.ICRU90 = 104.8,
    delta.R = dat$exact.delta,
    delta.ICRU90 = 2.470e-4,
    mu.st = dat$mu.st,
    L = dat$L)

  ############################
  # 10 keV
  ############################
  MeV <- 0.010 # Electron kinetic energy
  dat <- Sternheimer.delta.exact(MeV, dat)
  xx0 <- electronic.MSP.Bethe(MeV, dat, delta = 0)               # Compute MSP without density effect correction (delta = 0)
  xx  <- electronic.MSP.Bethe(MeV, dat, delta = dat$exact.delta) # Compute MSP with exact Sternheimer density correction

  df101 <- data.frame(
    MeV = MeV,
    I.eV = dat$I,
    rho=dat$exact.rho,
    MSP.R0 = xx0,
    MSP.R = xx,
    MSP.ICRU90 = 19.99,
    delta.R = dat$exact.delta,
    delta.ICRU90 = 2.634e-3,
    mu.st = dat$mu.st,
    L = dat$L)

  ############################
  # 100 keV
  ############################
  MeV <- 0.100 # Electron kinetic energy
  dat <- Sternheimer.delta.exact(MeV, dat)
  xx0 <- electronic.MSP.Bethe(MeV, dat, delta = 0)               # Compute MSP without density effect correction (delta = 0)
  xx  <- electronic.MSP.Bethe(MeV, dat, delta = dat$exact.delta) # Compute MSP with exact Sternheimer density correction

  df102 <- data.frame(
    MeV = MeV,
    I.eV = dat$I,
    rho=dat$exact.rho,
    MSP.R0 = xx0,
    MSP.R = xx,
    MSP.ICRU90 = 3.654,
    delta.R = dat$exact.delta,
    delta.ICRU90 = 4.047e-2,
    mu.st = dat$mu.st,
    L = dat$L)

  ############################
  # 800 keV
  ############################
  MeV <- 0.8 # Electron kinetic energy
  dat <- Sternheimer.delta.exact(MeV, dat)
  xx0 <- electronic.MSP.Bethe(MeV, dat, delta = 0)               # Compute MSP without density effect correction (delta = 0)
  xx  <- electronic.MSP.Bethe(MeV, dat, delta = dat$exact.delta) # Compute MSP with exact Sternheimer density correction

  df103 <- data.frame(
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


  ############################
  # 1 MeV
  ############################
  MeV <- 1.0 # Electron kinetic energy
  dat <- Sternheimer.delta.exact(MeV, dat)
  xx0 <- electronic.MSP.Bethe(MeV, dat, delta = 0)               # Compute MSP without density effect correction (delta = 0)
  xx  <- electronic.MSP.Bethe(MeV, dat, delta = dat$exact.delta) # Compute MSP with exact Sternheimer density correction

  df104 <- data.frame(
    MeV = MeV,
    I.eV = dat$I,
    rho=dat$exact.rho,
    MSP.R0 = xx0,
    MSP.R = xx,
    MSP.ICRU90 = 1.606,
    delta.R = dat$exact.delta,
    delta.ICRU90 = 0.7593,
    mu.st = dat$mu.st,
    L = dat$L)


  ############################
  # 2 MeV
  ############################
  MeV <- 2.0 # Electron kinetic energy
  dat <- Sternheimer.delta.exact(MeV, dat)
  xx0 <- electronic.MSP.Bethe(MeV, dat, delta = 0)               # Compute MSP without density effect correction (delta = 0)
  xx  <- electronic.MSP.Bethe(MeV, dat, delta = dat$exact.delta) # Compute MSP with exact Sternheimer density correction

  df105 <- data.frame(
    MeV = MeV,
    I.eV = dat$I,
    rho=dat$exact.rho,
    MSP.R0 = xx0,
    MSP.R = xx,
    MSP.ICRU90 = 1.586,
    delta.R = dat$exact.delta,
    delta.ICRU90 = 1.364,
    mu.st = dat$mu.st,
    L = dat$L)

  ############################
  # 10 MeV
  ############################

  MeV <- 10 # Electron kinetic energy
  dat <- Sternheimer.delta.exact(MeV, dat)
  xx0 <- electronic.MSP.Bethe(MeV, dat, delta = 0)               # Compute MSP without density effect correction (delta = 0)
  xx  <- electronic.MSP.Bethe(MeV, dat, delta = dat$exact.delta) # Compute MSP with exact Sternheimer density correction

  df106<- data.frame(
    MeV = MeV,
    I.eV = dat$I,
    rho=dat$exact.rho,
    MSP.R0 = xx0,
    MSP.R = xx,
    MSP.ICRU90 = 1.729,
    delta.R = dat$exact.delta,
    delta.ICRU90 = 3.381,
    mu.st = dat$mu.st,
    L = dat$L)

  ############################
  # 100 MeV
  ############################
  MeV <- 100 # Electron kinetic energy
  dat <- Sternheimer.delta.exact(MeV, dat)
  xx0 <- electronic.MSP.Bethe(MeV, dat, delta = 0)               # Compute MSP without density effect correction (delta = 0)
  xx  <- electronic.MSP.Bethe(MeV, dat, delta = dat$exact.delta) # Compute MSP with exact Sternheimer density correction

  df107 <- data.frame(
    MeV = MeV,
    I.eV = dat$I,
    rho=dat$exact.rho,
    MSP.R0 = xx0,
    MSP.R = xx,
    MSP.ICRU90 = 1.928,
    delta.R = dat$exact.delta,
    delta.ICRU90 = 7.624,
    mu.st = dat$mu.st,
    L = dat$L)



  ############################
  # 1000 MeV
  ############################
  MeV <- 1000 # Electron kinetic energy
  dat <- Sternheimer.delta.exact(MeV, dat)
  xx0 <- electronic.MSP.Bethe(MeV, dat, delta = 0)               # Compute MSP without density effect correction (delta = 0)
  xx  <- electronic.MSP.Bethe(MeV, dat, delta = dat$exact.delta) # Compute MSP with exact Sternheimer density correction

  df108 <- data.frame(
    MeV = MeV,
    I.eV = dat$I,
    rho=dat$exact.rho,
    MSP.R0 = xx0,
    MSP.R = xx,
    MSP.ICRU90 = 2.106,
    delta.R = dat$exact.delta,
    delta.ICRU90 = 12.22,
    mu.st = dat$mu.st,
    L = dat$L)


  df <- rbind(df100,df101,df102,df103,df104,df105,df106,df107,df108)

  ############################################################################
  # Main results
  ############################################################################

  # Excellent agreement with between the clanElectrons computations
  # and ICRU-90 values for water. MSP.R0 is the mass electronic
  # stopping power without density effect correction (delta=0).
  #
#  MeV I.eV   rho   MSP.R0    MSP.R MSP.ICRU90    sdelta.R delta.ICRU90
#  8e-01   81 2.265 1.694586 1.639646      1.640  0.6074777       0.6075
#  1e+00   81 2.265 1.671790 1.606030      1.606  0.7593199       0.7593
#  2e+00   81 2.265 1.694645 1.585502      1.586  1.3640823       1.3640
#  1e+01   81 2.265 1.989286 1.729347      1.729  3.3811044       3.3810
#  1e+02   81 2.265 2.512917 1.928154      1.928  7.6239905       7.6240
#  1e+03   81 2.265 3.042536 2.105602      2.106 12.2158180      12.2200

df
} # graphite (demo)


