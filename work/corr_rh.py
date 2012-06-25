from netCDF4 import *
from wrfpy import *
from mk_sounding_consts import consts
#--------------------
con = consts()
#-- funcs -----------
def interp(x, x1, x2, v1, v2):
  v = (v2-v1)/(x2-x1) * (x - x2) + v2
  return v
#--------------------
flag = 1
wrfinput = con.wrfexe_dir + "/wrfinput_d01"
n       = Dataset(wrfinput,"r","NETCDF")
lrh_nc  = wrf_rh(n, flag)[0,:]
lz_nc   = (n.variables["PHB"][0,:-1] + n.variables["PH"][0,:-1])/9.81
#q       = n.variables["QVAPOR"][0,:]
#theta   = n.variables["T"][0]  + 300.
#-- average ---------
lrh_nc  = mean(mean(lrh_nc, axis=2), axis=1)
lz_nc   = mean(mean(lz_nc , axis=2), axis=1)

#-- read input_sounding ---
#rhsound_org = "./input_sounding_rh2.txt"
#rhsound_mod = rhsound_org[:-4] + "_crr" + ".txt"
rhsound_org = con.sound_rh_org
rhsound_mod = con.sound_rh_crr

f = open(rhsound_org, "r")
lines    = f.readlines()
f.close()
#-- surface ------
ssurf        = lines[0]
ltemp        = map(float, ssurf.strip().split("\t"))
z_surf       = ltemp[0]
theta_surf   = ltemp[1]
rh_surf      = ltemp[2]
#-----------------
lines    = lines[1:]
nl_in    = len(lines)
#-- read rh org ---------
lz_sound      = []
ltheta_sound  = []
lrh_sound     = []
lu_sound      = []
lv_sound      = []
lrh_sound_mod = []
for i in range(nl_in):
  line   = map(float, lines[i].strip().split("\t"))
  lz_sound.append(line[0])
  ltheta_sound.append(line[1])
  lrh_sound.append(line[2])
  lu_sound.append(line[3])
  lv_sound.append(line[4])
#--- correct rh -----
n_nc      = len(lrh_nc)
n_sound    = nl_in
#-------------------
i_nc_start = 0
for i_sound in range(n_sound):
  rh_sound   =  lrh_sound[i_sound]
  z_sound    =  lz_sound[i_sound]
  for i_nc in range(i_nc_start, n_nc):
    if lz_nc[i_nc] >= z_sound:
      #---------
      if i_nc == 0:
        z_nc_pre    = z_surf
        rh_nc_pre   = rh_surf
      else:
        z_nc_pre    = i_nc -1
        rh_nc_pre   = lrh_nc[z_nc_pre]
        i_nc_start  = i_nc -1
      #---------
      z_nc_nxt      = lz_nc[i_nc]
      rh_nc_nxt     = lrh_nc[i_nc]
      break
    if ((i_nc == n_sound-1) & (lz_nc[i_nc] < z_sound)):
      #z_nc_pre      = lz_nc[i_nc-1]
      #z_nc_nxt      = lz_nc[i_nc]
      #rh_nc_pre     = lrh_nc[i_nc-1]
      #rh_nc_nxt     = lrh_nc[i_nc]
      z_nc_pre      = lz_nc[i_nc]
      z_nc_nxt      = z_sound
      rh_nc_pre     = lrh_nc[i_nc]
      rh_nc_nxt     = lrh_nc[i_nc]


  #-----------
  rh_nc     = interp(z_sound, z_nc_pre, z_nc_nxt, rh_nc_pre, rh_nc_nxt)
  drh       = rh_sound - rh_nc
  rh_sound  = rh_sound + drh
  lrh_sound_mod.append(rh_sound)

  print i_sound, drh, rh_sound, rh_nc
#---- write to file ---------------------------

#sout = "%s\t%s\t%s\t\t\n"%(lsurf[0], lsurf[1], qv_surf_in*100.)
sout = ssurf.strip() + "\n"
lout = []
for i in range(nl_in):
  zk     = lz_sound[i]
  theta  = ltheta_sound[i]
  rh     = lrh_sound_mod[i]     # [no dimension]
  u      = lu_sound[i]
  v      = lv_sound[i]
  stemp  = "%i\t%9.5f\t%13.10e\t%8.3f\t%8.3f\n"%(zk, theta, rh, u, v)
  sout   = sout + stemp
#
f = open(rhsound_mod, "w")
f.write(sout)
f.close()  



