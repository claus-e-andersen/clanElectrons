#
#
#
get.NIST.estar.data <- function(){

df.NIST.estar.SP <- read.table("C:\\data\\projects\\R\\8900-Theoretical-Dosimetry\\stopping-powers2\\NIST-estar-stopping-power-data.txt",header=TRUE,sep=";")

}

############################
get.NIST.estar <- function(MeV=c(0.1,1,10),what="MSP.el",id="water.liquid",data=df.NIST.estar.SP){
  df <- data
  ok <- df$id ==id
  df <- df[ok,]
  if(sum(ok)>0){
    # some data
    df %>%
      dplyr::select(MeV) ->
      xx
    xx <- xx[,1]

    df %>%
      dplyr::select(all_of(what)) ->
      yy

    yy <- yy[,1]
    zz <- approx(xx,yy,xout=MeV)$y
    zz
  }
}

######################
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
