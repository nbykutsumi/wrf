#-----------------------------------------------
# To: make a input_sounding data with
# regulated relative humidity (RH)
# from input_sounding_rhxxx.txt
#***********************************************
import os
from mk_sounding_consts import consts
#***********************************************
prodir       = "/home/utsumi/inst/wrf/work"
pro_mk_sound = prodir + "/mk_sounding.py"
pro_corr_rh  = prodir + "/corr_rh.py"
#--- consts ------------------------
con = consts()
wrfexe_dir   = con.wrfexe_dir
sounddir     = con.sounddir
sounding     = con.sounding
#--- step1 -------------------------
os.system("python %s 1"%(pro_mk_sound))

#--- step2 -------------------------
#
#
os.chdir(wrfexe_dir)
os.system("ln -sf %s ./input_sounding"%(sounding))
os.system("./ideal.exe")
os.system("cp -f ./wrfinput_d01 ./wrfinput_d01_org")
os.chdir(prodir)
#--- step3 -------------------------
#
# make a correction to input_sounding_rhxx.txt 
# here, input_sounding_rhxxx_crr.txt" is created
#---------
os.system("python %s"%(pro_corr_rh))

#--- step3--------------------------
#
# make a sounding data from corrected rh sounding
# here, input_sounding_xxx.txt is overwritten
#----------
os.system("python %s 2"%(pro_mk_sound))

