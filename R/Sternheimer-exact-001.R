#' @title Sternheimer.set.to.conductor
#' @description  Ascertain that material is treated as conductor in Sternheimer delta computations
#'
#'
#'  dat <- Sternheimer.set.to.conductor(dat)
#'
#' The binding energy of the outer subshell is set to zero
#' @export
Sternheimer.set.to.conductor <- function(dat){
# Set outer binding energy to zero
Evec <- dat$Evec.org
nlev <- length(Evec)
dat$Evec[nlev] <- 0
dat$fvec <- dat$fvec.org
dat$type <- "conductor"
dat
}

#' @title Sternheimer.set.to.insulation
#' @description  Ascertain that material is treated as insulator in Sternheimer delta computations
#'
#'  dat <- Sternheimer.set.to.insulation(dat)
#'
#' The binding energy of the outer subshell is NOT set to zero
#' @export
Sternheimer.set.to.insulator <- function(dat){
  # Set binding energies to original values
  dat$Evec <- dat$Evec.org
  dat$fvec <- dat$fvec.org
  dat$type <- "insulator"
  dat
}



#' @title Sternheimer.f.root.mu.st <-
#' @description  Helper function (mu.st) for exact Sternheimer density correction
#'
#' # Equation suitable for finding mu.st (eq. 4.29 in ICRU-90)
#' @export
Sternheimer.f.root.mu.st <- function(mu.st,dat){
# Helper function (mu.st) for exact Sternheimer density correction
# Created: Aug 7, 2022
# Revised: Aug 7, 2022
# Name   : Claus E. Andersen
# Modified version of eq. 4.29 in ICRU-90
# suitable for root finding: The requested
# x = mu.st will fullfill: f.root.mu.st(x) = 0
fvec <- dat$fvec
Evec <- dat$Evec
Ep <- dat$Ep
nlev <- dat$nlev
I <- dat$I
ans <- 0
fn <- fvec[nlev]
for(i in 1:(nlev)){
  if(fvec[i]>0){
  ans <- ans + fvec[i]* log( ( (mu.st*Evec[i])^2 + 2/3*fvec[i]*Ep^2)^0.5)
  }# if
}# loop
#ans <- ans + fn * log((fn)^0.5 * Ep)
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
# Revised: Aug 7, 2022
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

ans <- 0
for(i in 1:(nlev)){
  if(fvec[i]>0){
  ans <- ans + fvec[i] / ( (mu.st*Evec[i]/Ep)^2 + L^2 )
  }# if
}# loop

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
Sternheimer.delta.exact <- function(dat=NULL){
# Created: Aug 7, 2022
# Revised: Aug 7, 2022
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

if(is.null(dat)){
 dat <- list(
    plot.wanted = FALSE,
    MeV = 1000, # Kinetic energy of the electron
    nlev = 6,
    Z    = 13,
    A    = 26.98154,
    rho.density =  2.265,
    fvec = c(2/13, 2/13 ,2/13, 2/13, 2/13, 3/13),
    #Evec =c(  1564.0 , 121. , 77.0 , 77.0 , 10.62 , 5.986),
    Evec =c(  1564.0 , 121. , 77.0 , 77.0 , 10.62 ,      0), # Conductor
    I = 166.0)
  } # null dat

E0       <- 0.51099895000
dat$beta <- (1 - (E0/(E0+dat$MeV))^2)^0.5
dat$Ep   <- 28.8159 * (dat$rho.density *dat$Z/dat$A)^0.5
print(dat)


####################################################
# Part 1: Find mu.st (solve eq. 4.29 in ICRU-90)
####################################################
# Search for mu.st (coarse)
mu.st.min <- 0
mu.st.max <- 200
mu.st.N   <- 10000
mu.st.eps <- 0.00006

xx <- seq(mu.st.min,mu.st.max,length=mu.st.N)
yy <- Sternheimer.f.root.mu.st(xx,dat)

ok.pos <- yy>0
ok.neg <- yy<0

if(sum(ok.pos)> 0 & sum(ok.neg)>0){
# Search for mu.st (fine)
mu.st.min <- max(xx[ok.neg])
mu.st.max <- min(xx[ok.pos])
mu.st.N   <- 1000
mu.st.eps <- 0.0006

xx <- seq(mu.st.min,mu.st.max,length=mu.st.N)
yy <- Sternheimer.f.root.mu.st(xx,dat)

ok <- abs(yy) < mu.st.eps & xx>0
mu.st.root <- NA
if(sum(ok)>0){
  # Use regression to find final value (root)
  fm <- lm(xx[ok]~yy[ok])
  mu.st.root <- as.numeric(coefficients(fm)[1])
}

df <- data.frame(mu.st=xx,val=yy)

plt.mu.st <- lattice::xyplot(val ~ mu.st,
  data=df,
  panel=function(x,y,...){
    panel.xyplot(x,y,...)
    panel.abline(h=0,lty="dashed")
    panel.abline(v=mu.st.root,lty="dashed")
  }
)
if(dat$plot.wanted)(print(plt.mu.st))
}# both pos and neg

dat$mu.st <- mu.st.root

####################################################
# PART 2: Find L (solve eq. 4.28 in ICRU-90)
####################################################
# Search interval for L (coarse)
L.min <- -1
L.max <- 2000
L.N   <- 4000
L.eps <- 0.0006

xx <- seq(L.min,L.max,length=L.N)
yy <- Sternheimer.f.root.L(xx,dat)

ok.pos <- yy>0
ok.neg <- yy<0

if(sum(ok.pos)> 0 & sum(ok.neg)>0){
# Search interval for L (fine)
L.min <- max(xx[ok.pos])
L.max <- min(xx[ok.neg])
L.N   <- 10000
L.eps <- 0.0006

xx <- seq(L.min,L.max,length=L.N)
yy <- Sternheimer.f.root.L(xx,dat)

ok <- abs(yy) < L.eps & xx>0
L.root <- NA
if(sum(ok)>0){
  # Use regression to find final value (root)
  fm <- lm(xx[ok]~yy[ok])
  L.root <- as.numeric(coefficients(fm)[1])
}

df <- data.frame(L=xx,val=yy)
plt.L <- lattice::xyplot(val ~ L,
data=df,
panel=function(x,y,...){
  panel.xyplot(x,y,...)
  panel.abline(h=0,lty="dashed")
  panel.abline(v=L.root,lty="dashed")
}
)

if(dat$plot.wanted){print(plt.L)}
} # both pos and neg

dat$L  <- L.root

L     <- dat$L
mu.st <- dat$mu.st
Evec  <- dat$Evec
fvec  <- dat$fvec
beta  <- dat$beta
nlev  <- dat$nlev
Ep    <- dat$Ep

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
delta <- ans - L^2 * (1-beta^2)

# Done
dat$delta <- delta

if(dat$plot.wanted){clanLattice::print.trellis.plots(list(plt.mu.st,plt.L),1)}

dat
} #Sternheimer.delta.exact.function


#############################################################################


#' @title demo.Sternheimer.delta.exact
#' @export
demo.Sternheimer.delta.exact <- function(){
# Biniding energies can be found in
# Carlson (1975) Book: Photoelectron and Auger Spectroscopy
# p. 338. Table A1.A
# Binding Energies of Electrons in Free Atom (eV) : Z = 1-53
# See also: G4AtomicShells.cc (Geant4)

dat.Al <- list(
plot.wanted = FALSE,
MeV = 1000,
nlev = 6,
Z    = 13,
A    = 26.98154,
rho.density =  2.265,
fvec.org = c(2/13, 2/13 ,2/13, 2/13, 2/13, 3/13),
Evec.org = c( 1564.0 , 121.0, 77.0, 77.0, 10.62, 5.986),
I = 166.0)
dat.Al <- Sternheimer.set.to.conductor(dat.Al)

dat.C <- list(
plot.wanted = FALSE,
MeV = 1000,
nlev = 3,
Z    = 6,
A    = 12.011,
rho.density =  2.265,
fvec.org = c(2/6,2/6,2/6),
Evec.org = c(288.00, 16.59, 11.26),
I = 78)
dat.C <- Sternheimer.set.to.conductor(dat.C)


dat.Be <- list(
plot.wanted = !FALSE,
MeV = 1000,
nlev = 2,
Z    = 4,
A    = 9.01218,
rho.density =  1.848,
fvec.org = c(0.5,0.5),
Evec.org = c(115.0, 9.322),
I = 63.7)
dat.Be <- Sternheimer.set.to.insulator(dat.Be)

dat <- dat.Al
dat$MeV <- 10
dat <- Sternheimer.delta.exact(dat=dat)
dat
} # demo.Sternheimer.delta.exact

########################################################################################

#' @title demo.Sternheimer.delta.exact.plot
#' Delta vs. energy for two different models (Al: En=0 or not)
#' @export
demo.Sternheimer.delta.exact.plot <- function(){
# Created: July 29, 2022
# Revised: July 31, 2022
# Name:    Claus E. Andersen
print("This function should be run manually, line by line.")


dat.Al.model1 <- list(
plot.wanted = FALSE,
MeV = 1000,
nlev = 6,
Z    = 13,
A    = 26.98154,
rho.density =  2.265,
fvec.org = c(2/13, 2/13 ,2/13, 2/13, 2/13, 3/13),
Evec.org = c( 1564.0 , 121.0, 77.0, 77.0, 10.62, 5.986),
I = 166.0)

dat.Al.model1 <- Sternheimer.set.to.conductor(dat.Al.model1)
dat.Al.model2 <- Sternheimer.set.to.insulator(dat.Al.model1)

xx <- seq(-0.6,3,length=50)
ee <- 10^xx
yy1 <- 0 * xx
yy2 <- 0 * xx


dat <- dat.Al.model1
for(ii in 1:length(ee)){
  dat$MeV <- ee[ii]
  dat.out <- Sternheimer.delta.exact(dat=dat)
  yy1[ii] <- dat.out$delta
}

dat <- dat.Al.model2
for(ii in 1:length(ee)){
  dat$MeV <- ee[ii]
  dat.out <- Sternheimer.delta.exact(dat=dat)
  yy2[ii] <- dat.out$delta
}


df1 <- data.frame(MeV=ee, material = "Al (conductor model)", delta=yy1)
df2 <- data.frame(MeV=ee, material = "Al (insulator model)", delta=yy2)
df <- rbind(df1,df2)

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
} # demo.Sternheimer.delta.exact.plot
