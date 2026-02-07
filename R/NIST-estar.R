
#' @title df.NIST.estar.SP
#' @description  Extract data from the df.NIST.estar.SP data for given id (=name of element or material).
#' df.NIST.estar.SP is included in the package.
#' The data were downloaded from the NIST web site February 2-5, 2026.
#' # Each element has been downloaded from the NIST Estar data base and saves as an ASCII
#' file format like:

#'ESTAR: Stopping Powers and Range Tables for Electrons
#'
#'CALIFORNIUM                                                                 #
#'
#'Kinetic   Collision Radiative Total     CSDA      Radiation D. Effect
#'Energy    Stp. Pow. Stp. Pow. Stp. Pow. Range     Yield     Parameter
#'MeV       MeV cm2/g MeV cm2/g MeV cm2/g g/cm2
#'
#'1.000E-02 7.806E+00 2.175E-02 7.828E+00 9.638E-04 1.384E-03 1.044E-03
#'1.250E-02 6.848E+00 2.417E-02 6.872E+00 1.306E-03 1.737E-03 1.332E-03
#'1.500E-02 6.129E+00 2.621E-02 6.155E+00 1.691E-03 2.096E-03 1.631E-03
#'1.750E-02 5.568E+00 2.799E-02 5.596E+00 2.117E-03 2.458E-03 1.941E-03
#'2.000E-02 5.118E+00 2.956E-02 5.147E+00 2.584E-03 2.822E-03 2.263E-03
#'2.500E-02 4.437E+00 3.225E-02 4.469E+00 3.630E-03 3.554E-03 2.942E-03
#'3.000E-02 3.945E+00 3.453E-02 3.979E+00 4.818E-03 4.286E-03 3.669E-03
#'3.500E-02 3.571E+00 3.652E-02 3.608E+00 6.140E-03 5.017E-03 4.445E-03
#'4.000E-02 3.277E+00 3.829E-02 3.315E+00 7.587E-03 5.744E-03 5.270E-03
#'4.500E-02 3.039E+00 3.990E-02 3.079E+00 9.154E-03 6.468E-03 6.146E-03
#' ...
#'8.000E+02 1.554E+00 1.373E+02 1.389E+02 2.984E+01 9.504E-01 8.108E+00
#'9.000E+02 1.562E+00 1.547E+02 1.563E+02 3.052E+01 9.547E-01 8.336E+00
#'1.000E+03 1.568E+00 1.722E+02 1.737E+02 3.112E+01 9.583E-01 8.541E+00
#'
#' @format A data frame with 8240 rows and 15 variables:
#' \describe{
#'   \item{MeV}{Vector of energies (in MeV) for which the function should estimate the stopping power.}
#'   \item{MSP.el}{Electron (collision) mass stopping power (in MeV per g/cm2).}
#'   \item{MSP.rad}{Radiative mass stopping power (in MeV per g/cm2).}
#'   \item{MSP.tot}{Totalmass stopping power (in MeV per g/cm2).}
#'   \item{range.CSDA}{SSDA range (in g/cm2).}
#'   \item{rad.yield}{Radiation yield}
#'   \item{density.effect}{Density effect parameter}
#'   \item{type}{element (like "hydrogen") or material (like "water.liquid").}
#'   \item{id}{Name of element or material from NIST}
#'   \item{id2}{Name of elements from alternative source}
#'   \item{symbol}{Element symbols like "H", "He", "Br".}
#'   \item{rho}{Density (in g/cm3).}
#'   \item{I}{Mean excitation energy (eV)}
#'   \item{Z}{Element number (=0 for composite materials))}
#'   \item{source}{Source of data ("NIST.estar)"}
#'
#' }
#' @source Generated for package examples
"df.NIST.estar.SP"


#' @title df.NIST.xray.materials
#' @description NIST material data from X-Ray Mass Attenuation Coefficients database (not Estar)
#' # From https://physics.nist.gov/PhysRefData/XrayMassCoef/tab1.html
#' X-Ray Mass Attenuation Coefficients
#' Table 1. Material constants assumed in the present evaluations for elemental media.
#' Values are given for the ratio of atomic number-to-mass Z/A, the mean excitation energy I, and the density ρ.
#' Some density values are only nominal; those for Z = 85 and 87 were arbitrarily set to 10 in order to complete the calculations.
#' Variable 1: Z
#' Variable 2: Element
#' Variable 3: Name
#' Variable 4:  Z/A = Z.A.ratio
#' Variable 5: I	in eV
#' Variable 6: Density in g/cm3
#' Downloaded: February 7, 2026 by Claus E. Andersen
#'
#' @format A data frame with 8240 rows and 15 variables:
#' \describe{
#'   \item{Z}{Element number (=0 for composite materials))}
#'   \item{symbol}{Element symbols like "H", "He", "Br".}
#'   \item{id.xray}{Name of element or material from NIST}
#'   \item{Z.A.ratio}{Z/A}
#'   \item{I.xray}{Mean excitation energy (eV)}
#'   \item{rho.xray}{Density (in g/cm3).}
#' }
#' @source Generated for package examples
"df.NIST.xray.materials"


#
#
#' @title get.NIST.estar.data
#' @description  How the df.NIST.estar.SP data was created.
#'
#' @export
get.NIST.estar.data <- function(){
# This code is just to maintain a record of how the data was created.
# Do not run!
if(FALSE){ # FALSE BLOCK
  require(dplyr)
  require(dplyr)
  require(lattice)
  # Created: Feb. 2, 2026
  # Revised: Feb. 5, 2026
  # Revised: Feb. 7, 2026 I value for Z=93 + carbon (2 g/cm3) as a material (Z=0)
  # Revised: Feb. 7, 2026 Get Z/A + sybb as data)ols from NIST x-ray table (also included
  # Name: Claus E. Andersen

  # Each element has been downloaded from the NIST Estar data base and saves as an ASCII
  # file format like:

  #ESTAR: Stopping Powers and Range Tables for Electrons
  #
  #CALIFORNIUM                                                                 #
  #
  #Kinetic   Collision Radiative Total     CSDA      Radiation D. Effect
  #Energy    Stp. Pow. Stp. Pow. Stp. Pow. Range     Yield     Parameter
  #MeV       MeV cm2/g MeV cm2/g MeV cm2/g g/cm2
  #
  #1.000E-02 7.806E+00 2.175E-02 7.828E+00 9.638E-04 1.384E-03 1.044E-03
  #1.250E-02 6.848E+00 2.417E-02 6.872E+00 1.306E-03 1.737E-03 1.332E-03
  #1.500E-02 6.129E+00 2.621E-02 6.155E+00 1.691E-03 2.096E-03 1.631E-03
  #1.750E-02 5.568E+00 2.799E-02 5.596E+00 2.117E-03 2.458E-03 1.941E-03
  #2.000E-02 5.118E+00 2.956E-02 5.147E+00 2.584E-03 2.822E-03 2.263E-03
  #2.500E-02 4.437E+00 3.225E-02 4.469E+00 3.630E-03 3.554E-03 2.942E-03
  #3.000E-02 3.945E+00 3.453E-02 3.979E+00 4.818E-03 4.286E-03 3.669E-03
  #3.500E-02 3.571E+00 3.652E-02 3.608E+00 6.140E-03 5.017E-03 4.445E-03
  #4.000E-02 3.277E+00 3.829E-02 3.315E+00 7.587E-03 5.744E-03 5.270E-03
  #4.500E-02 3.039E+00 3.990E-02 3.079E+00 9.154E-03 6.468E-03 6.146E-03
  # ...
  #8.000E+02 1.554E+00 1.373E+02 1.389E+02 2.984E+01 9.504E-01 8.108E+00
  #9.000E+02 1.562E+00 1.547E+02 1.563E+02 3.052E+01 9.547E-01 8.336E+00
  #1.000E+03 1.568E+00 1.722E+02 1.737E+02 3.112E+01 9.583E-01 8.541E+00

  # We do not need this anymore
  #df.elements <- read.table("elementlist.txt",sep=",")
  #names(df.elements) <- c("Z","symbol","id2")
  #df.elements$id2 <- tolower(df.elements$id2)

  # Data fron NIST X-ray table (see header of file)
  df.xray.elements <- read.table("NIST-xray-materials.txt",skip=11)
  names(df.xray.elements) <- c("Z","symbol","id.xray","Z.A.ratio","I.xray","rho.xray")
  df.xray.elements$id2 <- tolower(df.xray.elements$id2)


  df.meta <- rbind(data.frame(id="air.dry"     ,fn="NIST-air-dry.txt",rho=1.20479E-03,I=85.7,Z=0),
                   data.frame(id="pmma"        ,fn="NIST-pmma.txt",rho=1.19,I=74.0,Z=0),
                   data.frame(id="alanine"     ,fn="NIST-alanine.txt",rho=1.41,I=71.9,Z=0),
                   data.frame(id="water.liquid",fn="NIST-water-liquid.txt",rho=1.0,I=75.0,Z=0),
                   data.frame(id="carbon",fn="NIST-carbon.txt",rho=2.0,I=81.0,Z=0),
                   #
                   data.frame(id="hydrogen",fn="NIST-hydrogen.txt",rho=8.37480E-05,I=19.2,Z=1),
                   data.frame(id="helium",fn="NIST-helium.txt",rho=1.66322E-04,I=41.8,Z=2),
                   data.frame(id="lithium",fn="NIST-lithium.txt",rho=5.34000E-01,I=40.0,Z=3),
                   data.frame(id="beryllium",fn="NIST-beryllium.txt",rho=1.848,I=63.7,Z=4),
                   data.frame(id="boron",fn="NIST-boron.txt",rho=2.37,I=76,Z=5),
                   data.frame(id="graphite",fn="NIST-graphite.txt",rho=1.7,I=78.0,Z=6),
                   data.frame(id="nitrogen",fn="NIST-nitrogen.txt",rho=1.16528E-03,I=82.0,Z=7),
                   data.frame(id="oxygen",fn="NIST-oxygen.txt",rho=1.33151E-03,I=95.0,Z=8),
                   data.frame(id="fluorine",fn="NIST-fluorine.txt",rho=1.58029E-03,I=115,Z=9),
                   data.frame(id="neon",fn="NIST-neon.txt",rho=8.38505E-04,I=137,Z=10),
                   data.frame(id="sodium",fn="NIST-sodium.txt",rho=9.71000E-01,I=149,Z=11),
                   data.frame(id="magnesium",fn="NIST-magnesium.txt",rho=1.74000E+00,I=156,Z=12),
                   data.frame(id="aluminum",fn="NIST-aluminum.txt",rho=2.69890,I=166.0,Z=13),
                   data.frame(id="silicon",fn="NIST-silicon.txt",rho=2.33,I=173.0,Z=14),
                   data.frame(id="phosphorus",fn="NIST-phosphorus.txt",rho=2.2,I=173.0,Z=15),
                   data.frame(id="sulfur",fn="NIST-sulfur.txt",rho=2.0,I=180,Z=16),
                   data.frame(id="chlorine",fn="NIST-chlorine.txt",rho=2.99473E-03,I=174.0,Z=17),
                   data.frame(id="argon",fn="NIST-argon.txt",rho=1.66201E-03,I=188.0,Z=18),
                   data.frame(id="potassium",fn="NIST-potassium.txt",rho=8.62000E-01,I=190,Z=19),
                   data.frame(id="calcium",fn="NIST-calcium.txt",rho=1.55,I=191.0,Z=20),
                   data.frame(id="scandium",fn="NIST-scandium.txt",rho=2.989,I=216,Z=21),
                   data.frame(id="titanium",fn="NIST-titanium.txt",rho=4.54,I=233,Z=22),
                   data.frame(id="vanadium",fn="NIST-vanadium.txt",rho=6.11,I=245.0,Z=23),
                   data.frame(id="chromium",fn="NIST-chromium.txt",rho=7.18,I=257.0,Z=24),
                   data.frame(id="manganese",fn="NIST-manganese.txt",rho=7.44,I=272.0,Z=25),
                   data.frame(id="iron",fn="NIST-iron.txt",rho=7.874,I=286.0,Z=26),
                   data.frame(id="cobalt",fn="NIST-cobalt.txt",rho=8.9,I=297.0,Z=27),
                   data.frame(id="nickel",fn="NIST-nickel.txt",rho=8.902,I=311,Z=28),
                   data.frame(id="copper",fn="NIST-copper.txt",rho=8.96,I=322.0,Z=29),
                   data.frame(id="zinc",fn="NIST-zinc.txt",rho=7.133,I=330,Z=30),
                   data.frame(id="gallium",fn="NIST-gallium.txt",rho=5.904,I=334.0,Z=31),
                   data.frame(id="germanium",fn="NIST-germanium.txt",rho=5.323,I=350.0,Z=32),
                   data.frame(id="arsenic",fn="NIST-arsenic.txt",rho=5.73,I=347.0,Z=33),
                   data.frame(id="selenium",fn="NIST-selenium.txt",rho=4.50,I=348.0,Z=34),
                   data.frame(id="bromine",fn="NIST-bromine.txt",rho=7.07218E-03,I=343.0,Z=35),
                   data.frame(id="krypton",fn="NIST-krypton.txt",rho=3.47832E-03,I=352.0,Z=36),
                   data.frame(id="rubidium",fn="NIST-rubidium.txt",rho=1.5320,I=363,Z=37),
                   data.frame(id="strontium",fn="NIST-strontium.txt",rho=2.54,I=366.0,Z=38),
                   data.frame(id="yttrium",fn="NIST-yttrium.txt",rho=4.469,I=379.0,Z=39),
                   data.frame(id="zirconium",fn="NIST-zirconium.txt",rho=6.506,I=393.0,Z=40),
                   data.frame(id="niobium",fn="NIST-niobium.txt",rho=8.57,I=417.0,Z=41),
                   data.frame(id="molybdenum",fn="NIST-molybdenum.txt",rho=1.02200E+01,I=424.0,Z=42),
                   data.frame(id="technetium",fn="NIST-technetium.txt",rho=1.15000E+01,I=428.0,Z=43),
                   data.frame(id="ruthenium",fn="NIST-ruthenium.txt",rho=1.24100E+01,I=441.0,Z=44),
                   data.frame(id="rhodium",fn="NIST-rhodium.txt",rho=1.24100E+01,I=449.0,Z=45),
                   data.frame(id="palladium",fn="NIST-palladium.txt",rho=1.20200E+01,I=470.0,Z=46),
                   data.frame(id="silver",fn="NIST-silver.txt",rho=10.5,I=470.0,Z=47),
                   data.frame(id="cadmium",fn="NIST-cadmium.txt",rho=8.65,I=469.0,Z=48),
                   data.frame(id="indium",fn="NIST-indium.txt",rho=7.31,I=488.0,Z=49),
                   data.frame(id="tin",fn="NIST-tin.txt",rho=7.31,I=488.0,Z=50),
                   data.frame(id="antimony",fn="NIST-antimony.txt",rho=6.691,I=487.0,Z=51),
                   data.frame(id="tellurium",fn="NIST-tellurium.txt",rho=6.24,I=485.0,Z=52),
                   data.frame(id="iodine",fn="NIST-iodine.txt",rho=4.93000,I=491,Z=53),
                   data.frame(id="xenon",fn="NIST-xenon.txt",rho=5.48536E-03,I=482,Z=54),
                   data.frame(id="cesium",fn="NIST-cesium.txt",rho=1.8730,I=488,Z=55),
                   data.frame(id="barium",fn="NIST-barium.txt",rho=3.5,I=491,Z=56),
                   data.frame(id="lanthanum",fn="NIST-lanthanum.txt",rho=6.154,I=501,Z=57),
                   data.frame(id="cerium",fn="NIST-cerium.txt",rho=6.657,I=523,Z=58),
                   data.frame(id="praseodymium",fn="NIST-praseodymium.txt",rho=6.71,I=535,Z=59),
                   data.frame(id="neodymium",fn="NIST-neodymium.txt",rho=6.9,I=546,Z=60),
                   data.frame(id="promethium",fn="NIST-promethium.txt",rho=7.22,I=560,Z=61),
                   data.frame(id="samarium",fn="NIST-samarium.txt",rho=7.46,I=574,Z=62),
                   data.frame(id="europium",fn="NIST-europium.txt",rho=5.243,I=580,Z=63),
                   data.frame(id="gadolinium",fn="NIST-gadolinium.txt",rho=7.90040,I=591,Z=64),
                   data.frame(id="terbium",fn="NIST-terbium.txt",rho=8.229,I=614,Z=65),
                   data.frame(id="dysprosium",fn="NIST-dysprosium.txt",rho=8.55,I=628,Z=66),
                   data.frame(id="holmium",fn="NIST-holmium.txt",rho=8.795,I=650,Z=67),
                   data.frame(id="erbium", fn="NIST-erbium.txt", rho=9.066,I=658,Z=68),
                   data.frame(id="thulium",fn="NIST-thulium.txt",rho=9.321,I=674,Z=69),
                   data.frame(id="ytterbium",fn="NIST-ytterbium.txt",rho=6.73,I=684,Z=70),
                   data.frame(id="lutetium",fn="NIST-lutetium.txt",rho=9.84,I=694,Z=71),
                   data.frame(id="hafnium",fn="NIST-hafnium.txt",rho=1.331e1,I=705,Z=72),
                   data.frame(id="tantalum",fn="NIST-tantalum.txt",rho=1.66540E+01,I=718,Z=73),
                   data.frame(id="tungsten",fn="NIST-tungsten.txt",rho=1.93000E+01,I=727,Z=74),
                   data.frame(id="rhenium",fn="NIST-rhenium.txt",rho=2.102e1,I=736,Z=75),
                   data.frame(id="osmium",fn="NIST-osmium.txt",rho=2.25700E+01,I=746,Z=76),
                   data.frame(id="iridium",fn="NIST-iridium.txt",rho=2.24200E+01,I=757,Z=77),
                   data.frame(id="platinum",fn="NIST-platinum.txt",rho=2.14500E+01,I=790,Z=78),
                   data.frame(id="gold",fn="NIST-gold.txt",rho=1.93200E+01,I=790.0,Z=79),
                   data.frame(id="mercury",fn="NIST-mercury.txt",rho=1.35460E+01,I=800.0,Z=80),
                   data.frame(id="thallium",fn="NIST-thallium.txt",rho=1.17200E+01,I=810.0,Z=81),
                   data.frame(id="lead",fn="NIST-lead.txt",rho=1.13500E+01,I=823.0,Z=82),
                   data.frame(id="bismuth",fn="NIST-bismuth.txt",rho=9.747,I=823.0,Z=83),
                   data.frame(id="polonium",fn="NIST-polonium.txt",rho=9.32,I=830.0,Z=84),
                   data.frame(id="astatine",fn="NIST-astatine.txt",rho=9.32,I=825.0,Z=85),
                   data.frame(id="radon",fn="NIST-radon.txt",rho=9.06618E-03,I=794.0,Z=86),
                   data.frame(id="francium",fn="NIST-francium.txt",rho=1,I=827.0,Z=87),
                   data.frame(id="radium",fn="NIST-radium.txt",rho=5,I=826.0,Z=88),
                   data.frame(id="actinium",fn="NIST-actinium.txt",rho=1.007E+01,I=841.0,Z=89),
                   data.frame(id="thorium",fn="NIST-thorium.txt",rho=1.17200E+01,I=847.0,Z=90),
                   data.frame(id="protactinium",fn="NIST-protactinium.txt",rho=	1.53700E+01,I=878.0,Z=91),
                   data.frame(id="uranium",fn="NIST-uranium.txt",rho=1.89500E+01,I=890.0,Z=92),
                   data.frame(id="neptunium",fn="NIST-neptunium.txt",rho=2.02500E+01,I=902.0,Z=93),
                   data.frame(id="plutonium",fn="NIST-plutonium.txt",rho=1.98400E+01,I=921.0,Z=94),
                   data.frame(id="americium",fn="NIST-americium.txt",rho=1.36700E+01,I=934.0,Z=95),
                   data.frame(id="curium",fn="NIST-curium.txt",rho=1.35100E+01,I=939.0,Z=96),
                   data.frame(id="berkelium",fn="NIST-berkelium.txt",rho=1.40000E+01,I=952.0,Z=97),
                   data.frame(id="californium",fn="NIST-californium.txt",rho=10,I=966.0,Z=98)
  )

  # The I-value in estar (794 eV) for Z=86 (radon) seems to be wrong.
  # We therefore will use the NIST xray data (see below).
  # However, we keep rho and I from the NIST Estar data base as
  # rho.estar and I.estar, respectively.
  names(df.meta) <- c("id","fn","rho.estar","I.estar","Z")

  #df.meta %>%
  #left_join(df.elements) ->
  #df.meta


  df.meta %>%
    left_join(df.xray.elements) ->
    df.meta



  df <- NULL
  pn <- ".\\data"

  for(i in 1:nrow(df.meta)){
    print(i)
    fn <- paste(pn,"\\",df.meta$fn[i],sep="")
    print(fn)
    df0 <- read.table(fn,skip=9)
    names(df0) <- c("MeV","MSP.el","MSP.rad","MSP.tot","range.CSDA","rad.yield","density.effect")
    df0 %>%
      mutate(type = ifelse(df.meta$Z[i]==0,"material","element")) %>%
      mutate(id = df.meta$id[i]) %>%
      mutate(id2 = df.meta$id.xray[i]) %>%
      mutate(symbol = df.meta$symbol[i]) %>%
      mutate(rho = df.meta$rho.estar[i]) %>%
      mutate(I = df.meta$I.estar[i]) %>%
      mutate(Z = df.meta$Z[i]) %>%
      mutate(Z.A.ratio = df.meta$Z.A.ratio[i]) %>%
      mutate(source="NIST.estar") ->
      df0

    if(is.null(df)){df <- df0} else {df <- rbind(df,df0)}
  } # meta loop


  df %>%
    filter(!id==id2)

  df %>%
    arrange(Z,id,I) -> df

  head(df)

  tail(df)

  head(df %>% filter(Z %in% 1:92))

  df %>%
    filter(!is.na(Z.A.ratio)) %>%
    tail(.)

  #filter(!is.na(I.estar)) %>%
  #filter(Z %in% 1:98)%>%
  #filter(abs(rho.estar - rho) > 0.01) %>%
  #filter(MeV==1)


  write.table(df,"export100.txt",row.names=FALSE,sep=";")
  write.table(df.xray.elements,"xray-export200.txt",row.names=FALSE,sep=";")

  #require(lattice)
  #xyplot(log(MSP.el) ~ log(MeV),groups=rho<0.01,data=df,type="p")

  str(df)


  # Or read the finished file
  df.NIST.estar.SP <- read.table("C:\\data\\projects\\R\\8900-Theoretical-Dosimetry\\stopping-powers2\\NIST-estar-stopping-power-data.txt",header=TRUE,sep=";")
  df.NIST.xray.materials <- read.table("C:\\data\\projects\\R\\8900-Theoretical-Dosimetry\\stopping-powers2\\NIST-xray-materials-data.txt",header=TRUE,sep=";")

  # usethis::use_data(df.NIST.estar.SP,overwrite=TRUE)
  # usethis::use_data(df.NIST.xray.materials,overwrite=TRUE)
} # FALSE BLOCK
  } # how to produce the data


#' @title get.estar.from.id
#' @description  Extract data from the df.NIST.estar.SP data for given id (=name of element or material).
#' df.NIST.estar.SP is included in the package.
#'
#' @param MeV = vector of kinetic energies of the electron in MeV
#' @param what ="MSP.el", MSP.rad", "MSP.tot", "CSDA", "rad.yield" or "density.effect"
#' @param id = "hydrogen", "helium" .. "californium", "water.liquid", "air.dry", "pmma", "alanine", "graphite"
#' @param data = date frame with data (df.NIST.estar.SP)
#' @param rule = rule for approx: 1=no extrapolation, 2=use nearest value for extrapolation
#' @details Notes:
#' The  MSP's are in units of MeV pr. g/cm2.
#' The I is in eV
#' The density (rho) is in g/cm3
#' @export
get.estar.from.id <- function(MeV=c(0.1,1,10),what="MSP.el",id="water.liquid",data=df.NIST.estar.SP,rule=1){
  # Estimate NIST Estar data at given MeV energies using the approx interpolation function.
  # This version select the data based on given id
  if(length(unique(what))>1){print("Problem (get.estar.from.id): Multiple what's were defined. Only first one will be selected.")}
  if(length(unique(id))>1){print("Problem (get.estar.from.id): Multiple is's were defined. Only first one will be selected.")}
  df <- data
  id <- toupper(id)
  ok <- toupper(df$id) %in% id[1]
  df <- df[ok,]
  if(sum(ok)>0){
    # some data
    df |>
      dplyr::select(MeV) ->
      xx

    xx <- xx[,1]

    df |>
      dplyr::select(all_of(what)) ->
      yy

    yy <- yy[,1]
    zz <- approx(xx,yy,xout=MeV,rule=rule)$y
    zz
  }
} # end get.estar.from.id


#' @title get.estar.from.Z
#' @description  Extract data from the df.NIST.estar.SP data for given Z.
#' df.NIST.estar.SP is included in the package.
#'
#' @param MeV = vector of kinetic energies of the electron in MeV
#' @param what ="MSP.el", MSP.rad", "MSP.tot", "CSDA", "rad.yield" or "density.effect"
#' @param Z = element number (1,2,3 .. 98) (composite materials have Z=0)
#' @param data = date frame with data (df.NIST.estar.SP)
#' @param rule = rule for approx: 1=no extrapolation, 2=use nearest value for extrapolation
#' @details Notes:
#' The  MSP's are in units of MeV pr. g/cm2.
#' The I is in eV
#' The density (rho) is in g/cm3
#' @export
get.estar.from.Z<- function(MeV=c(0.1,1,10),what="MSP.el",Z=1,data=df.NIST.estar.SP,rule=1){
  # Estimate NIST Estar data at given MeV energies using the approx interpolation function.
  # This version select the data based on given Z
  if(length(unique(what))>1){print("Problem (get.estar.from.id): Multiple what's were defined. Only first one will be selected.")}
  if(length(unique(Z))>1){print("Problem (get.estar.from.Z): Multiple Z's were defined. Only first one will be selected.")}
  df <- data
  ok <- df$Z %in% Z[1]
  df <- df[ok,]
  if(sum(ok)>0){
    # some data
    df |>
      dplyr::select(MeV) ->
      xx

    xx <- xx[,1]

    df |>
      dplyr::select(all_of(what)) ->
      yy

    yy <- yy[,1]
    zz <- approx(xx,yy,xout=MeV,rule=rule)$y
    zz
  }
} # end get.estar.from.Z

########################################


######################
#' @title NIST.estar.demo
#' @description  Demonstration  of how to use get.estar.from.Z and get.estar.from.id
#'
#' @export
NIST.estar.demo <- function(){

  get.estar.from.id(id="lead",what=c("MSP.rad"))

  ff <- get.estar.from.id(MeV=seq(0.2,2,length=100),what="MSP.el",id="air.dry")
  plot(seq(0.2,2,length=100),ff)

  plt <- lattice::xyplot(log(MSP.el) ~ log(MeV) | paste("Z=",Z),data=df.NIST.estar.SP,type="l")
  print(plt)

  print(tail(df.NIST.estar.SP))

  print('get.estar.from.Z(MeV=c(0.1,1,10,20), Z=82, what="MSP.rad") : ')
    tmp <- get.estar.from.Z(MeV=c(0.1,1,10,20), Z=82, what="MSP.rad")
  print(tmp)

  print('get.estar.from.id(MeV=c(0.1,1,10,20), id="water.liquid", what="MSP.el") : ')
  tmp <- get.estar.from.id(MeV=c(0.1,1,10,20), id="water.liquid", what="MSP.el")
  print(tmp)

  print('get.estar.from.id(MeV=c(0.1,1,10,20), id="air.dry", what="MSP.el") : ')
  tmp <- get.estar.from.id(MeV=c(0.1,1,10,20), id="air.dry", what="MSP.el")
  print(tmp)

} # demo end



