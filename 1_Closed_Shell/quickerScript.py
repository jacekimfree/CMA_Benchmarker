import shutil as sh
import os
import sys

here = os.getcwd()
dirs = os.listdir(here)

for i in dirs:
    if os.path.isdir(os.getcwd() + "/" + i):
        os.chdir(i)
        os.chdir("CCSD_T_TZ")
        # os.chdir("Disps_B3LYP_6-31G_2df,p_")
        os.chdir("Disps_CCSD_T_DZ")
        files = os.listdir(os.getcwd())
        for j in range(len(files)):
            if os.path.exists(os.getcwd()+"/output."+str(j+1)+".dat"):
                print(j)
                os.remove("output."+str(j+1)+".dat")
        # for j in range(len(files)):
            # if os.path.isdir(i):
                # print(j)
                # os.remove("output."+str(i+1)+".dat")
        os.chdir('..')
        os.chdir('..')
        os.chdir('..')
