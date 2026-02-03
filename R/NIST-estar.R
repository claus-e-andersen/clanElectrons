#
#
#
#' @title get.NIST.estar.data
#' @description  How the df.NIST.estar.SP data was created.
#'
#' @export
get.NIST.estar.data <- function(){

df.NIST.estar.SP <- read.table("C:\\data\\projects\\R\\8900-Theoretical-Dosimetry\\stopping-powers2\\NIST-estar-stopping-power-data.txt",header=TRUE,sep=";")

}
#' @title get.NIST.estar
#' @description  Extract data from the df.NIST.estar.SP data frame included
#' with the package.
#'
#' @param MeV = vector of kinetic energies of the electron in MeV
#' @param what ="MSP.el", MSP.rad", "MSP.tot", "CSDA", "rad.yield" or "density.effect"
#' @param id = "water.liquid", "air.dry", "pmma", "alanine"
#' @param data = date frame with data (df.NIST.estar.SP)
#' @details Notes:
#' The  MSP's are in units of MeV pr. g/cm2.
#' The I is in eV
#' The density (rho) is in g/cm3
#' @export
############################
get.NIST.estar <- function(MeV=c(0.1,1,10),what="MSP.el",id="water.liquid",data=df.NIST.estar.SP){
  df <- data
  ok <- df$id ==id
  df <- df[ok,]
  if(sum(ok)>0){
    # some data
    df |>
      dplyr::select(MeV) ->
      xx
    xx <- xx[,1]

    df |>
      dplyr::select(dplyr::all_of(what)) ->
      yy

    yy <- yy[,1]
    zz <- stats::approx(xx,yy,xout=MeV)$y
    zz
  }
}

######################
#' @title get.NIST.estar.demo
#' @description  Demonstration  of how to use get.NIST.estar
#'
#' @export
get.NIST.estar.demo <- function(){

ff <- get.NIST.estar(id="pmma")
print(ff)
ff <- get.NIST.estar(id="alanine")
print(ff)
ff <- get.NIST.estar(id="air.dry")
print(ff)


ff <- get.NIST.estar(MeV=seq(0,2,length=100),what="density.effect",id="water.liquid")
print(ff)
plot(ff)

}

########################
