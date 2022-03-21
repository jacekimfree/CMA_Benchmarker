import shutil as sh
import os
import sys

here = os.getcwd()
dirs = os.listdir(here)

for i in dirs:
    wd = here + "/" + i
    if os.path.isdir(wd):
        os.chdir(i+"/CCSD_T_TZ/Disps_B3LYP_6-31G_2df,p_/")
        files = os.listdir(os.getcwd())
        if os.path.exists(os.getcwd()+"/old/"):
            sh.rmtree("/old")
        if not "/fc_int.dat" in files:
            print(os.getcwd() + "needs fc_int.dat")
        for j in range(len(files)):
            if os.path.exists(os.getcwd()+"/output."+str(j+1)+".dat"):
                os.remove("output."+str(j+1)+".dat")
        os.chdir(here)
