from sounding import *
esflag = 1
lsound = sounding.sounding_rh2q_zlev(esflag)

p_surf   = lsound[0]
th_surf  = lsound[1]
qv_surf  = lsound[2]
rh_surf  = lsound[3]
azk      = lsound[4]
ap       = lsound[5]
atheta   = lsound[6]
arho     = lsound[7]
au       = lsound[8]
av       = lsound[9]
aqv      = lsound[10]
arh      = lsound[11]
atk      = lsound[12]
nl_in    = lsound[13]

#-------------------------
# pressure, Pa --> hPa
p_surf   = p_surf *0.01
ap       = ap*0.01
#-------------------------
# mixing ratio,  kg/kg --> g/kg
qv_surf  = qv_surf * 1000.
aqv      = aqv * 1000.
#-------------------------

soname   = "./input_sounding_myrun2.txt"
csvname  = "./check_input_sounding_myrun2.csv"
#-----------------------------
sout = "%i\t%6.2f\t%7.4f\n"%(p_surf, th_surf, qv_surf)
for i in range(nl_in):
  stemp  = "%i\t%6.2f\t%14.9e\t%.3f\t%.3f\n"%(azk[i], atheta[i], aqv[i], au[i], av[i])
  sout   = sout + stemp
#
f = open(soname, "w")
f.write(sout)
f.close()

#--- csv file for check ------
sout  = ""
ltemp = [p_surf, th_surf, qv_surf]
sout  = ",".join(map(str,ltemp)) + "\n"
sout  = sout + "zk, theta, qv, pm, tk, RH, rho" + "\n"
for i in range(nl_in):
  ltemp    = [azk[i], atheta[i], aqv[i], ap[i], atk[i], arh[i], arho[i]]
  sout     = sout + ",".join(map(str, ltemp)) + "\n"



f = open(csvname, "w")
f.write(sout)
f.close()






