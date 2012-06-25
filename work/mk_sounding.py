from sep_module_initialize_ideal import *
from mk_sounding_consts import consts
import sys
#-----------------------------------------------
con = consts()
#-----------------------------------------------
def cal_es(t, flag):
  #------
  # calc es (Pa)
  #------
  if flag == 1:
    a      = 7.5
    b      = 237.3
  elif flag == -1:
    a      = 9.5
    b      = 265.5
  #---
  tdeg   = t - 273.16
  es = 6.11 * 10**(a * tdeg / (b + tdeg)) * 100.0
  return es
#----------------
def cal_qs(p, t, flag):
  #  qs (kg kg-1)
  #  p  (Pa)
  #  t  (K)
  epsi   = 0.62185
  es     = cal_es(t, flag)
  qs     = epsi * es / (p - es)
  return qs
#-----------------------------------------------
stepflag     = int(sys.argv[1])
if stepflag  == 1:
  #soundfile   = "./input_sounding_rh.txt"
  #soundfile   = "./input_sounding_rh2.txt"
  soundfile    = con.sound_rh_org
elif stepflag == 2:
  #soundfile   = "./input_sounding_rh_crr.txt"
  #soundfile   = "./input_sounding_rh2_crr.txt"
  soundfile    = con.sound_rh_crr
#-----------------------------------------------
#soname   = "./input_sounding_myrun2.txt"
#csvname  = "./check_input_sounding_myrun2.csv"
soname   = con.sounding
csvname  = con.check_csv

esflag   = 1
rd       = 287.0
rd_cp    = 0.286
g        = 9.8
#- read sounding -------------------------------
f           = open(soundfile, "r")
lsound      = f.readlines()
f.close()
#
n           = len(lsound)
nl_in       = n -1
lsurf       = map(float, lsound[0].strip().split("\t"))
p_surf_in   = lsurf[0]   # (hPa)
th_surf_in  = lsurf[1]   # (K)
rh_surf_in  = lsurf[2]   # (%)
#-- convert units ---
p_surf_in   = p_surf_in * 100.0   #(Pa)
#--------------------
t_surf_in  = th_surf_in * (100000.0 / p_surf_in) ** rd_cp
qv_surf_in = cal_qs( p_surf_in, t_surf_in, esflag) * rh_surf_in  #[kg kg-1]
tv_surf_in = t_surf_in * (1 + 0.61*qv_surf_in)
#--------------------

#
lzk         = []
ltheta      = []
lrh         = []
lu          = []
lv          = []
#
for i in range(1,n):
  ltemp     = map(float, lsound[i].strip().split("\t"))
  lzk.append( ltemp[0])
  ltheta.append( ltemp[1])
  lrh.append( ltemp[2])
  lu.append( ltemp[3])
  lv.append( ltemp[4]) 
#-- calc qv[0] ---------------------
dz   =  lzk[0] 
dp   =  p_surf_in * g / (rd * tv_surf_in) *dz
p0   =  p_surf_in - dp
qv0  = cal_qs(p0, t_surf_in, esflag) * lrh[0] # (kg/kg)
#-----------------------------------------------
lqv_temp      = [qv0]
for nl_temp in range(2, nl_in+1):
#for nl_temp in range(2,5+1):
  lqv_temp    = lqv_temp + [9999.]
  #
  nl_in       =  nl_temp
  lzk_temp    = lzk[:nl_temp]
  ltheta_temp = ltheta[:nl_temp]
  lrh_temp    = lrh[:nl_temp]
  rh_target   = lrh_temp[-1]
  #------------------------

  #sep_module_initialize_ideal.init_domain_rk(p_surf_in, th_surf_in, qv_surf_in, lzk, ltheta_temp, lqv_temp, nl_max, nl_in)

  qv_top  =  sep_module_initialize_ideal.mk_lqv(p_surf_in, th_surf_in, rh_surf_in, lzk_temp, ltheta_temp, lqv_temp, rh_target)
  print "inpy, qv_top, rh",qv_top, lrh_temp[-1]
  lqv_temp[-1] = qv_top
lqv  = lqv_temp
print "lqv",lqv
#---- write to file ---------------------------

sout = "%s\t%s\t%s\t\t\n"%(lsurf[0], lsurf[1], qv_surf_in*1000.)
lout = []
for i in range(nl_in):
  zk     = lzk[i]
  theta  = ltheta[i]
  qv     = lqv[i] * 1000.  # [g/kg]
  u      = lu[i]
  v      = lv[i]
  #stemp  = "%i\t%9.5f\t%16.11f\t%8.3f\t%8.3f\n"%(zk, theta, qv, u, v)
  #stemp  = "%i\t%9.5f\t%11.6f\t%8.3f\t%8.3f\n"%(zk, theta, qv, u, v)
  #stemp  = "%i\t%9.5f\t%11.7f\t%8.3f\t%8.3f\n"%(zk, theta, qv, u, v)
  #stemp  = "%i\t%9.5f\t%11.8f\t%8.3f\t%8.3f\n"%(zk, theta, qv, u, v)
  #stemp  = "%i\t%9.5f\t%13.9f\t%8.3f\t%8.3f\n"%(zk, theta, qv, u, v)
  stemp  = "%i\t%9.5f\t%13.10e\t%8.3f\t%8.3f\n"%(zk, theta, qv, u, v)
  #stemp  = "%i\t%9.5f\t%11.6e\t%8.3f\t%8.3f\n"%(zk, theta, qv, u, v)
  sout   = sout + stemp
#
f = open(soname, "w")
f.write(sout)
f.close()

#--- csv file for check ------
sout  = ""
ltemp = [lsurf[0], lsurf[1], qv_surf_in*1000.]
sout  = ",".join(map(str,ltemp)) + "\n"
sout  = sout + "zk, theta, qv" + "\n"
for i in range(nl_in):
  zk     = lzk[i]
  theta  = ltheta[i]
  qv     = lqv[i] * 100.  # [g/kg]
  ltemp    = [zk, theta, qv]
  sout     = sout + ",".join(map(str, ltemp)) + "\n"

f = open(csvname, "w")
f.write(sout)
f.close()

