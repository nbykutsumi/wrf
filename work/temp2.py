from netCDF4 import *
n = Dataset("./wrfout_d01_0001-01-01_00:00:00.myrun2.nc","r","NETCDF")
# calculate absolute temperature (K)
# Rd / CP   = 0.28571 (dimensionless)
# Pair      = P + PB
# theta     = Tk * (100000./Pair)^(Rd/Cp)
# --> Tk    = theta * (Pair / 100000.)^(R\/Cp)
#------------------------------------
#Rd_Cp       = 0.28571
Rd          = 287.0
Cp          = 7.0*Rd / 2.0
Rd_Cp       = Rd / Cp
Pair        = n.variables["P"][:] + n.variables["PB"][:]
theta       = n.variables["T"][:] + 300.
tk          = theta * (( Pair/100000. )**(Rd_Cp))

# calculate Relative humidity (%) with respect to liquid water
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
# q       = QVAPOR
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
Tair      = tk
es        = 10.0*SVP1*exp(SVP2* (Tair - SVPT0)/(Tair - SVP3))
QVS       = EP_3 * es / (0.01*Pair - (1.0 - EP_3)*es)

rh        = ma.masked_greater(QV/QVS, 1.0).filled(1.0)
rh        = ma.masked_less(rh, 0.0).filled(0.0)
rh        = 100.0 * rh




#epsi       = 0.622
#q          = n.variables["QVAPOR"][:]
#Pair       = n.variables["P"][:] + n.variables["PB"][:]
#theta      = n.variables["T"][:] + 300.
#e_sat      = 0.611 * exp( 17.2694 * (tk - 273.16) / (tk - 35.86) ) * 1000.0  # [Pa]
#e          = q * Pair / (q + epsi)
##
#rh         = e / e_sat *100.



## calculate Relative humidity (%) with respect to liquid water
##------------------------------------
## input variables 
## - T: perturbationj potential temperature (theta-t0)";
## - QVAPOR: Water vapor mixing ratio (kg kg-1)
## - P: perturbation pressure
## - PB: base state pressure
##------------------------------------
## epsi    = 0.622
## Rd / CP = 0.28571 (dimensionless)
## Pair    = P + PB
## theta   = T + 300
## Tair    = theta * ( Pair /100000.)**(Rd/Cp)
## e_sat   = 0.611 * exp( 17.2694 * (Tair - 273.16) / (Tair - 35.86) )  # Teten's formula
## q       = QVAPOR
## e       = q * Pair / (q + epsi)
##------------------------------------
## rh      = e / e_sat * 100.
##------------------------------------
#epsi       = 0.622
#Rd_Cp      = 0.28571
#q          = n.variables["QVAPOR"][:]
#Pair       = n.variables["P"][:] + n.variables["PB"][:]
#theta      = n.variables["T"][:] + 300.
#Tair       = theta * ( Pair/100000. )**(Rd_Cp)
#e_sat      = 0.611 * exp( 17.2694 * (Tair - 273.16) / (Tair - 35.86) )
#e          = q * Pair / (q + epsi)
##
#rh         = e / e_sat *100.
#return rh
