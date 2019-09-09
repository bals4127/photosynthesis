#' @title Calculate photosynthetic discrimination for C and O isotopes.
#'
#'@description Calculate photosynthetic discrimination for C and O isotopes using 4 cycle caliberations and one licor connection.
#' One can define upto two different combinations of 4 cycle caliberations.
#'
#' @param dat clean data file for TDl (.csv).
#' @param site_line_seq_1 site seqence used in TDL. This function expect 7 site sequnece: 4 (20,21,22,23)- caliberations, 1(3)- caliberation gas,  2 (7,8)- Licor reference and sample
#' @param site_line_seq_2 secod site seqence (if any) used in TDL. This function expect 7 site sequnece: 4 (20,21,22,24)- caliberations, 1(3)- caliberation gas,  2 (7,8)- Licor reference and sample
#' @param licor_ref name licor reference line number
#' @param licor_samp name licor sample line number
#' @param L3_tot_CO2 Total CO2 conecntration of caliberation gas for line 3
#' @param L3_d13_pdb d13C (pdb)of CO2 in caliberation gas for line 3
#' @param d18O_L3_raw d18O (pdb) of CO2 conecntration of caliberation gas for line 3
#' @param L3_d18O_SMOW d18O (vsmow) of CO2 conecntration of caliberation gas for line 3
#' @param d13pdb_mix d13C (pdb) of CO2 conecntration of mixing gas
#' @param d18Ovsmow_mixing_tank d18O (vsmow) of CO2 conecntration of mixing gas
#' @param f_other 13C fractionation asscoiated with other isotopologus species
#' @param f_other_o 18O fractionation asscoiated with other isotopologus species
#' @param Rvpdb_C the 13C/12C ratio of the reference standard, PDB
#' @param Rvpdb_O the 18O/16O ratio of the reference standard, PDB
#' @param R18Osmow the 18O/16O ratio of the reference standard, SMOW
#' @param . ddply function parameater
#' @return Lists with raw calcuations, clean output for discrimination and caliberations used in the calculations.
#' @importFrom plyr ddply
#' @importFrom stats reshape
#' @export
cal_tdl<- function(dat, site_line_seq_1= c(20,21,22,23,3,7,8),
                    site_line_seq_2= c(20,21,23,24,3,7,8),licor_ref= 7, licor_samp= 8,
                    L3_tot_CO2= 365.89, L3_d13_pdb = -37.68,
                    d18O_L3_raw= -20.5, L3_d18O_SMOW = 20.126552536249999,
                    d13pdb_mix = -40.88, d18Ovsmow_mixing_tank = 5.2230191347249963,
                    f_other = 0.00474, Rvpdb_C = 0.0111797,
                    f_other_o = 0.1185, Rvpdb_O = 0.01185, R18Osmow =0.0020052){
  dat<- dat[!duplicated(dat), ] ##remove duplicated
  RSCPDB3<- (L3_d13_pdb/1000+1)*Rvpdb_C+1 # (Eq 1)
  L3_tot_CO2_correct<- L3_tot_CO2*(1-f_other) #  (Eq 2)
  concA_avg3<- L3_tot_CO2_correct/(RSCPDB3)  # (Eq 3)
  concB_avg3<- L3_tot_CO2_correct - concA_avg3  # (Eq 4)
  measured_C16O2_L3 <- L3_tot_CO2 * (1- Rvpdb_O)/(1+2*R18Osmow*(1+L3_d18O_SMOW/1000)) # (Eq 6)Conc of 16O in CO2
  measured_C18O16O_L3 <- L3_tot_CO2*(1-Rvpdb_O)- measured_C16O2_L3   #   (Eq 7) Conc of 18O in CO2
  cal_const<- data.frame(RSCPDB3,L3_tot_CO2_correct,concA_avg3,
                         concB_avg3, measured_C16O2_L3,  measured_C18O16O_L3) ## Cal_constant
  dat$raw_tot_CO2 <- with(dat, Conc12CO2_Avg+ Conc13CO2_Avg)/(1-f_other)
  ##calculate raw delta
  dat$raw_delta <- with(dat,  1000*(((Conc13CO2_Avg/Conc12CO2_Avg)/Rvpdb_C) -1))
  assign_id <- function(dat, site_seq= site_line_seq_1){
    ###########################################################
    dat<- subset(dat, RECORD > 0) #Cleans of first junk cycle
    #################################################################
    #########     The following line may need to be changed based on your ramps used ###
    tdl.start <- which.min(dat$SiteOutput == site_seq[1]) #Starts at a relevant reference site.
    #########     The following line may need to be changed based on which site the LICOR is on ###
    tdl.end <- max(which(dat$SiteOutput == site_seq[7])) #Ends the data frame at the end of a cycle.
    #Begin and end with the good stuff
    dat<- dat[1:tdl.end,]
    ##start id with one
    id<-1
    #add empty column in dataframe
    dat$grpid<-0
    ##put in loop to add ID for repeating sites
    for (i in 1:nrow(dat)){
      print(i)
      if ((dat$SiteOutput[i]==site_seq[1] & !is.na(dat$SiteOutput[i])) && (dat$SiteOutput[i+1]==site_seq[2] & !is.na(dat$SiteOutput[i+1])) &&
          (dat$SiteOutput[i+2]==site_seq[3] & !is.na(dat$SiteOutput[i+2])) && (dat$SiteOutput[i+3]==site_seq[4] & !is.na(dat$SiteOutput[i+3]))&&
          (dat$SiteOutput[i+4]==site_seq[5] & !is.na(dat$SiteOutput[i+4])) &&(dat$SiteOutput[i+5]==site_seq[6]& !is.na(dat$SiteOutput[i+5])) &&
          (dat$SiteOutput[i+6]==site_seq[7] & !is.na(dat$SiteOutput[1+6])) )
      {
        dat$grpid[i]<-id
        dat$grpid[i+1]<-id
        dat$grpid[i+2]<-id
        dat$grpid[i+3]<-id
        dat$grpid[i+4]<-id
        dat$grpid[i+5]<-id
        dat$grpid[i+6]<-id
        id <- id+1
      }
    }
    print(c(i,id))
    return(dat)
  }
  ## Subzero function for each cycle
  sub_zero<- function (cycle){
    ## Substract site 20 from other sites
    sub_zero_C12<- with(cycle, (Conc12CO2_Avg)- (Conc12CO2_Avg[1])) # Column G
    sub_zero_C13<- with(cycle, (Conc13CO2_Avg)- (Conc13CO2_Avg[1])) # Column H
    sub_zero_tot_CO2<- (sub_zero_C12 + sub_zero_C13)/(1 - f_other) # Column I
    ## 13C: 13C16O16O  12C : 12C16O16O and 18O : 12C18O16O
    C16O2_raw<-  with(cycle, (raw_tot_CO2  + Conc18O_Avg)/(1-f_other_o))  # Column Y
    C18O2_raw <- cycle$Conc18O_Avg  # Column X
    subtract_zero_O16<- (C16O2_raw)- (C16O2_raw[1]) # Column AA
    subtract_zero_O18<- (C18O2_raw)- (C18O2_raw[1]) # Column AB
    subtract_zero_total_for_O2<- (subtract_zero_O16 + subtract_zero_O18)/(1-f_other_o) # Column AB
    sub_zero_dat<- data.frame(cycle, sub_zero_C12, sub_zero_C13, sub_zero_tot_CO2,
                              C16O2_raw,  C18O2_raw , subtract_zero_O16,subtract_zero_O18,subtract_zero_total_for_O2)
    return(sub_zero_dat)
  }
  fun_gain<- function(cycle){
    if (nrow(cycle)<7){return(cycle)} else{
      gain12<- concA_avg3/cycle[cycle$SiteOutput==3, "sub_zero_C12" ] # Column J
      gain13<- concB_avg3/cycle[cycle$SiteOutput==3, "sub_zero_C13" ] # Column K
      gain_12CO2<- with(cycle, sub_zero_C12* gain12)  # Column L
      czt_13CO2<-   gain_12CO2 *((d13pdb_mix/1000)+1) * Rvpdb_C # Column M, note that here it is using delta from mixing tank
      gain_16_O <- measured_C16O2_L3/cycle[cycle$SiteOutput==3, "subtract_zero_O16" ] # Column AD
      gain_18_O <- measured_C18O16O_L3/cycle[cycle$SiteOutput==3, "subtract_zero_O18" ] # Column AE
      czt_C16O <- with(cycle, subtract_zero_O16 * gain_16_O)  # Column AF
      czt_C18O16O <- czt_C16O *((d18Ovsmow_mixing_tank * 0.0020052/1000)+0.0020052)*2 # Column AG
      gain_dat<- data.frame(cycle, gain12, gain13, gain_12CO2, czt_13CO2, gain_16_O,gain_18_O,czt_C16O, czt_C18O16O )
      return(gain_dat)}}
  corr_cal_c<-  function(cycle){
    moddat<-cycle[1:4,]  ## Subset first four sites, mixing sites
    quad<- with(moddat, lm(czt_13CO2 ~ sub_zero_C13 +I(sub_zero_C13 ^2))) #fit quadratic model
    options(digits = 22)  ## this will collect model parameaters to 22 decimal points for more precision
    #Collect output of model fit in small dataframe
    a = as.numeric(summary(quad)$coefficients[1,1]) # Column N
    b = as.numeric(summary(quad)$coefficients[2,1]) # Column O
    c = as.numeric(summary(quad)$coefficients[3,1]) # Column P
    fit<- data.frame(a,b,c)
    #Correct values fro 12 and 13 CO2
    corr13<- fit[, "a"] + fit[, "b"]*cycle$sub_zero_C13  +fit[, "c"]*cycle$sub_zero_C13^2 # Column Q
    corr12<- cycle$gain_12CO2  # Column R , Yes it is correct values includes gain +raw CO2
    tmf_13CO2 <- (corr13+corr12)/ (1-f_other)   # Column T
    d13<- (((corr13/corr12)/Rvpdb_C)-1) * 1000  # Column U Note that unit is in per mil

    quad_O<- with(moddat, lm(czt_C18O16O ~ subtract_zero_O18+I(subtract_zero_O18^2)))
    #quad<- with(moddat, lm(subtract_zero_O18 ~ czt_C18O16O+I(czt_C18O16O^2)))  ##Bala Edits function corrected
    fit_O<- data.frame(a = summary(quad_O)$coefficients[1,1],
                       b = summary(quad_O)$coefficients[2,1],
                       c = summary(quad_O)$coefficients[3,1]) ## Column AH, AI and AJ
    corr18<- fit_O[, "a"] + fit_O[, "b"]*cycle$subtract_zero_O18  +  fit_O[, "c"]*cycle$subtract_zero_O18^2 #Column AK
    corr16<- cycle$czt_C16O # Column AL
    tmf_C18OO <- (corr18+corr16)/ (1-0.01185)   # Column AN
    d18<- (0.5*((corr18/corr16)/0.0020052)-1)*1000
    corr_C<- data.frame(cycle,corr12,corr13, tmf_13CO2, d13,  corr18,  corr16, tmf_C18OO, d18)
    return(corr_C)}
  grouprecs<-function(df){
    ref<- as.numeric(licor_ref)
    samp<- as.numeric(licor_samp)
    df$SiteOutput[df$SiteOutput==ref] <- 31
    df$SiteOutput[df$SiteOutput==samp] <- 32
    ## for universal way we need to replace those site by ideal one so that
    #we can keep following equations common after reshaping
    id<-1
    df$grpid<-0
    #df$tmf<-(df$corr12+df$corr13)/(1-0.00474)
    for (i in 1:nrow(df)-1){

      if ((df$SiteOutput[i]==31) && (df$SiteOutput[i+1]==32)){
        df$grpid[i]<-id
        df$grpid[i+1]<-id
        id<-id+1
      }
    }
    dfab<-subset(df,grpid>0)
    dfab$SiteOutput[dfab$SiteOutput==31] <- "Lic_ref"
    dfab$SiteOutput[dfab$SiteOutput==32] <- "Lic_samp"

    dfab1<-reshape(dfab[,c("DateTime","grpid","SiteOutput","corr12","corr13","d13", "d18","tmf_13CO2","tmf_C18OO")],
                   v.names=c("corr12","corr13","d13", "d18", "tmf_13CO2","tmf_C18OO"  ),
                   idvar="grpid",ids="DateTime",timevar="SiteOutput",direction="wide")
    # dfab$p<-dfab$tmf.19/(dfab$tmf.19-dfab$tmf.20)
    dfab1$p<-dfab1$tmf_13CO2.Lic_ref/(dfab1$tmf_13CO2.Lic_ref-dfab1$tmf_13CO2.Lic_samp)
    dfab1$p18O<-dfab1$tmf_C18OO.Lic_ref/(dfab1$tmf_C18OO.Lic_ref-dfab1$tmf_C18OO.Lic_samp)
    #dfab1$D13<-(dfab1$p*(((dfab1$d13.32/1000)-(dfab1$d13.31/1000))/ (1+(dfab1$d13.32/1000)-(dfab1$p*((dfab1$d13.32/1000)-(dfab1$d13.31/1000))))))*1000
    #dfab1$D18_O <- (dfab1$p *( (dfab1$d18.32/1000) - (dfab1$d18.31/1000))) / (1+ (dfab1$d18.32/1000) - (dfab1$p *  ( (dfab1$d18.32/1000) - (dfab1$d18.31/1000))) ) *1000
    ## same equations
    dfab1$D13<-((dfab1$p*(dfab1$d13.Lic_samp-dfab1$d13.Lic_ref))/ (1000+ dfab1$d13.Lic_samp-(dfab1$p*(dfab1$d13.Lic_samp-dfab1$d13.Lic_ref))))*1000
    dfab1$D18_O <- ((dfab1$p18O *(dfab1$d18.Lic_samp -dfab1$d18.Lic_ref)) / (1000+ dfab1$d18.Lic_samp - (dfab1$p18O *(dfab1$d18.Lic_samp - dfab1$d18.Lic_ref))))*1000
    return(dfab1)
  }
  assn_id <- assign_id(dat, site_seq= site_line_seq_1)
  non_assn<- droplevels(subset(assn_id, grpid== 0 ))
  if(nrow(non_assn)>10){
    warning("You have different site sequence,make sure site_line_seq_1 and site_line_seq_2 are correct,
            if you are sure ingore this message!!!")
    assn_2 <- assign_id(dat= non_assn, site_seq= site_line_seq_2)
    assn_2 <- assn_2[assn_2$grpid > 0, ]
    assn_1 <- droplevels(subset(assn_id, grpid > 0 ))

    zero_cal_1<- ddply(assn_1, "grpid", fun=sub_zero)
    gain_cal_1<- ddply(zero_cal_1, "grpid",fun=fun_gain)
    correct_c_1<- ddply(gain_cal_1, "grpid", fun=corr_cal_c)
    isotope_1<- grouprecs(correct_c_1)

    zero_cal_2<- ddply(assn_2, "grpid", fun=sub_zero)
    gain_cal_2<- ddply(zero_cal_2, "grpid",fun=fun_gain)
    correct_c_2<- ddply(gain_cal_2,"grpid",fun=corr_cal_c)
    isotope_2<- grouprecs(correct_c_2)
    isotope<- rbind(isotope_1, isotope_2)
    correct_c<- rbind(correct_c_1, correct_c_2)}
  else{
    assn <- droplevels(subset(assn_id, grpid > 0 ))
    zero_cal<- ddply(assn, "grpid", fun=sub_zero)
    gain_cal<- ddply(zero_cal, "grpid",fun=fun_gain)
    correct_c<- ddply(gain_cal, "grpid", fun=corr_cal_c)
    isotope<- grouprecs(correct_c)}
  l <- list()
  l$isotope<- isotope [order(isotope$DateTime),]
  l$raw_cal<- correct_c[order(correct_c$DateTime),]
  l$calibrations<- cal_const
  class(l) <- "TDL"
  return(l)
}
