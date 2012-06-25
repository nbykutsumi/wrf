from netCDF4 import *
from wrfpy import *

flag = 1
n = Dataset("./wrfinput_d01","r","NETCDF")
rh = wrf_rh(n, flag)[0,:]
q= n.variables["QVAPOR"][0,:]
z= (n.variables["PHB"][0,:-1] + n.variables["PH"][0,:-1])/9.81
theta = n.variables["T"][0]  + 300.
