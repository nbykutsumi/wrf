from numpy import *
def wrf_tk(n):
  # calculate absolute temperature (K)
  # Rd / CP   = 0.28571 (dimensionless)
  # Pair      = P + PB
  # theta     = Tk * (100000./Pair)^(Rd/Cp)
  # --> Tk    = theta * (Pair / 100000.)^(R\/Cp)
  #------------------------------------
  Rd          = 287.0
  Cp          = 7.0*Rd / 2.0
  Rd_Cp       = Rd / Cp
  Pair        = n.variables["P"][:] + n.variables["PB"][:]
  theta       = n.variables["T"][:] + 300.
  tk          = theta * (( Pair/100000. )**(Rd_Cp))
  return tk
#********************************************************************** 
def wrf_rh(n, flag):
  # calculate Relative humidity (%) with respect to liquid water
  # this function use "Tetens formula"
  #------------------------------------
  # input variables 
  # - T: perturbationj potential temperature (theta-t0)";
  # - QVAPOR: Water vapor mixing ratio (kg kg-1)
  # - P: perturbation pressure
  # - PB: base state pressure
  # # es is calculated with respect to : flag=1 water, flag=-1 ice, flag=0 mix
  #------------------------------------
  # epsi    = 0.622
  # Rd / CP = 0.28571 (dimensionless)
  # Pair    = P + PB   # (Pa)
  # theta   = T + 300  # (K)
  # Tair    = theta * ( Pair /100000.)**(Rd/Cp)    # (K)
  # e_sat   = 0.01 * 6.1078 * 10^( a * Tdeg / (b + Tdeg) )  # Teten's formula
  # q       = QVAPOR   # (kg kg-1), water vapor mixing ratio
  # e       = q * Pair / (q + epsi)
  #------------------------------------
  # rh      = e / e_sat * 100.
  #------------------------------------
  EP_3      = 0.622
  QV        = n.variables["QVAPOR"][:]
  Pair      = n.variables["P"][:] + n.variables["PB"][:]
  Tair      = wrf_tk(n)
  #-- calc es -----------------------------
  a_liq = 7.5
  b_liq = 237.3
  a_ice = 9.5
  b_ice = 265.3
  rTliq = 273.15  #   0 deg.C
  rTice = 250.15  # -23 deg.C
  #
  Tdeg     =  Tair- 273.16
  #  
  if (flag == 1):
    es = 100. * 6.1078 * 10.0**(a_liq * Tdeg/ (b_liq + Tdeg))
  elif (flag == -1):
    es = 100. * 6.1078 * 10.0**(a_ice * Tdeg/ (b_ice + Tdeg))
  elif (flag == 0):
    if ( Tdeg >= Tdegliq):
      es = 100. * 6.1078 * 10.0**(a_liq * Tdeg/ (b_liq + Tdeg))
    elif ( Tdeg <= Tdegice ):
      es = 100. * 6.1078 * 10.0**(a_ice * Tdeg/ (b_ice + Tdeg))
    else:
      es_liq = 100.* 6.1078 * 10.0**(a_liq * Tdeg/ (b_liq + Tdeg))
      es_ice = 100.* 6.1078 * 10.0**(a_ice * Tdeg/ (b_ice + Tdeg)) 
      es = ((Tdeg - Tdegice)*es_liq + (Tdegliq - Tdeg)*es_ice)/(Tdegliq - Tdegice)
  #----------------------------------------
  QVS       = EP_3 * es / (Pair - es)     # [kg kg-1]
  rh        = ma.masked_greater(QV/QVS, 1.0).filled(1.0)
  rh        = ma.masked_less(rh, 0.0).filled(0.0)
  rh        = rh
  return rh

#********************************************************************** 
def wrf_rh_ncl(n):
  # calculate Relative humidity (%) with respect to liquid water
  # this function follows the ncl function "wrf_rh"
  #------------------------------------
  # input variables 
  # - T: perturbationj potential temperature (theta-t0)";
  # - QVAPOR: Water vapor mixing ratio (kg kg-1)
  # - P: perturbation pressure
  # - PB: base state pressure
  #------------------------------------
  # epsi    = 0.622
  # Rd / CP = 0.28571 (dimensionless)
  # Pair    = P + PB
  # theta   = T + 300
  # Tair    = theta * ( Pair /100000.)**(Rd/Cp)
  # e_sat   = 0.611 * exp( 17.2694 * (Tair - 273.16) / (Tair - 35.86) )  # Teten's formula
  # q       = QVAPOR   # [kg kg-1]
  # e       = q * Pair / (q + epsi)
  #------------------------------------
  # rh      = e / e_sat * 100.
  #------------------------------------
  SVP1      = 0.6112
  SVP2      = 17.67
  SVP3      = 29.65
  SVPT0     = 273.15
  R_D       = 287.0
  R_V       = 461.6
  EP_2      = R_D/R_V
  EP_3      = 0.622
  QV        = n.variables["QVAPOR"][:]
  Pair      = n.variables["P"][:] + n.variables["PB"][:]
  Tair      = wrf_tk(n)
  es        = 10.0*SVP1*exp(SVP2* (Tair - SVPT0)/(Tair - SVP3))
  QVS       = EP_3 * es / (0.01*Pair - (1.0 - EP_3)*es)
  
  rh        = ma.masked_greater(QV/QVS, 1.0).filled(1.0)
  rh        = ma.masked_less(rh, 0.0).filled(0.0)
  rh        = rh
  return rh
