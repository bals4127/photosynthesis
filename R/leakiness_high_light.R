#' @title Leakiness calculations in C4 photosynthesis under high light measurments.
#'
#' @description Calculates leakiness under high light assuming that the,
#' 1. bundle-sheath CO2 concentration much higher than the CO2 concentration in mesophyll cells,
#' 2. High rate CO2 hydration i.e. Vp/Vh = 0 and Vc = A+Rd,
#' 3. No photorespiration Vo= 0, and
#' 4. Fractionation during day respiration = 0.
#'
#' Eqn referes to Ubierna et al., 2018 (See references).
#
#' @param data dataframe with variables
#' @param varnames List of names of variables in the dataset (see Details).
#' @param e CO2 fractionation during day respiration, assumed  = 0‰
#' @param a Weighter 12C/13 fractionation for diffusion across the boundary layer and stomata in series (4.4‰).
#' @param b3 12C/13C fractionation during Rubisco carboxylation (29‰).
#' @param s 12C/13C fractionation during CO2 diffusion in the boundary layer (2.9‰)
#' @param ab 12C/13C fractionation during CO2 leackage out bundle-sheath cell assuming there is no HCO3 leackage (1.8‰).
#' @param am Summed 12C/13C fractionation during liquid-phase diffusion and dissolution of CO2 (1.8‰)
#' @param ... Further arguments (ingore at the moment).
#' @return A dataframe with components:
#' \item{t}{ternary correction factor (‰); Eqn 9}
#' \item{b3_bar}{(‰); Eqn 62}
#' \item{b4_bar}{(‰); Eqn 63}
#' \item{e_prime}{12C/13C fractionation during decarboxylation including the effect of
#' a respiratoy substrate isotipically distinct from recent photosynthate (‰); Eqn 28}
#' \item{leakiness}{unitless; Eqn 52 and 64}
#'
#' @details List of variables and their units need to suply in \code{varnames}.
#' \code{unique_id}{: This is a unique id for each row. Makes easy to combine output with other dataframes.}
#' \code{Anet}{: Net rate of CO2 assimilation; umol m^-2s^-1,}
#' \code{Rd}{: Rate of respiration in light; umol m^-2s^-1,}
#' \code{CibyC}{: Ratio of intercelluar CO2 concentration to the ambient,}
#' \code{gm}{: Measoplyll cell CO2 conductance' umol m^-2s^-1 Pa^-1,}
#' \code{D13}{: Observed 13C photsynthetic discrimination (Eqn5); ‰,}
#' \code{Ci_Pa}{: CO2 partial pressure inside the leaf; Pa,}
#' \code{Ca_Pa}{: CO2 partial pressure in the ambient air; Pa,}
#' \code{CL_Pa}{: CO2 partial pressure at the leaf surface; Pa,}
#' \code{E}{: Transpiration rate; mol m^-2s^-1,}
#' \code{Cond_CO2}{: Total conductance to diffusion of CO2 in air; mol m^-2 s^-1,}
#' \code{Tleaf}{: Leaf temperature; degree C,}
#' \code{d13_growth_air}{: delta13C of the CO2 at the plant growth environment,}
#' \code{d13_measureair}{: delta13C of the CO2 in the air used during gas-exchange measurement (Licor reference line).}
#'
#' @references
#' Ubierna N, Holloway-Phillips M-M, Farquhar GD (2018) Using Stable Carbon Isotopes to
#' Study C3 and C4 Photosynthesis: Models and Calculations. In S Covshoff, ed, Photosynthesis: Methods and Protocols. Springer New York, New York, NY, pp 155–196
#'
#' von Caemmerer S, Ghannoum O, Pengelly JJL, Cousins AB (2014) Carbon isotope discrimination as a tool to explore C4 photosynthesis. J Exp Bot 65: 3459–3470.
#'
#' Farquhar G (1983) On the nature of carbon isotope discrimination in C4 species. Functional Plant Biol 10: 205–226
#' @export
leakiness_hl<- function(data, varnames= c(unique_id= "unique_id",Anet= "Photo", Rd = "Rd", CibyCa= "Ci.Ca", gm= "gm",
                                      D13= "D13", Ci_Pa= "Ci_Pa", Ca_Pa= "Ca_Pa", CL_Pa= "CL_Pa",
                                      E= "E", Cond_CO2 = "CndCO2", Tleaf= "Tleaf" ,d13_growth_air= "dgrowth_air",
                                      d13_measure_air = "d13.31"), e= 0,a=4.4, b3=29, s=1.8, ab=2.9,am= 1.8,...){## mKe new dataframe with varible names to avoid error due to different user names
  unique_id<- data[, varnames[["unique_id"]]]
  Photo<- data[, varnames[["Anet"]]]
  Rd<- data[, varnames[["Rd"]]]
  Ci.Ca<- data[, varnames[["CibyCa"]]]
  gm<- data[, varnames[["gm"]]]
  D13<- data[, varnames[["D13"]]]
  Ci_Pa<- data[, varnames[["Ci_Pa"]]]
  Ca_Pa<- data[, varnames[["Ca_Pa"]]]
  CL_Pa<- data[, varnames[["CL_Pa"]]]
  E<- data[, varnames[["E"]]]
  Cond_CO2<- data[, varnames[["Cond_CO2"]]]
  Tleaf<- data[, varnames[["Tleaf"]]]
  d13_growth_air<- data[, varnames[["d13_growth_air"]]]
  d13_measure_air<- data[, varnames[["d13_measure_air"]]]
  #______________________________________________________________________
  eprime  <- e + d13_measure_air - d13_growth_air
  ## Calculate a'
  ap_n <- ((1+(ab/1000))*(Ca_Pa- CL_Pa)) + ((1+a/1000)*(CL_Pa - Ci_Pa))
  ap_d<-  (Ca_Pa- Ci_Pa)
  ap<- ((ap_n/ap_d)-1)*1000
  b4_prime<-  (-(9.483*1000)/ (273 + Tleaf) ) + 23.89 +2.2
  b4_bar<-  b4_prime - (eprime*0.5* Rd)/(Photo+ 0.5* Rd)
  #calculate b3p
  b3_bar<- b3 - eprime*((Rd/(Photo+Rd)) - (0.5* Rd/(Photo + 0.5 * Rd)))
  ##caculate ternary term
  t<-(1+ap/1000)* E / (2 * Cond_CO2)
  by_t<- 1/(1+t)
  cal_t<- (1-t)/(1+t)
  ### Calculate leakines using high light assumptions
  leak<-  (cal_t*D13-by_t*ap- (am- b4_bar)*(Photo/(Ca_Pa * gm))-(b4_bar- ap/(1+t))*Ci.Ca )/((b3_bar-s)*(Ci.Ca- Photo/(Ca_Pa*gm)))
  out_data<- data.frame(unique_id= unique_id,e_prime= eprime,a_prime= ap,b4_bar= b4_bar, b3_bar= b3_bar, t= t,leakiness=leak)
  return(out_data)
  }
