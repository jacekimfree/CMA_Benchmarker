import shutil as sh
import os
import sys

here = os.getcwd()
dirs = os.listdir(here)

for i in dirs:
    if os.path.isdir(os.getcwd() + "/" + i):
        # print(i.split("_"))
        j = i.split("_")[0]
        os.chdir(i)
        os.chdir("CCSD_T_TZ")
        os.chdir("Disps_CCSD_T_DZ")
        print(i)
        print(j)
        # for k in os.listdir(os.getcwd()):
            # sh.rmtree(k)
        # if os.path.exists('/home/vulcan/nlk36701/1_projects/5_MH/2_G2_Test_Set/CMA1/1_CMA1/' + j + '/DispsInit'):
            # for k in os.listdir('/home/vulcan/nlk36701/1_projects/5_MH/2_G2_Test_Set/CMA1/1_CMA1/' + j + '/DispsInit'):
                # if os.path.isdir('/home/vulcan/nlk36701/1_projects/5_MH/2_G2_Test_Set/CMA1/1_CMA1/' + j + '/DispsInit/' + k) and int(j) < 102 and not os.path.exists(os.getcwd()+'/'+k):
                    # print(k)
                    # sh.copytree('/home/vulcan/nlk36701/1_projects/5_MH/2_G2_Test_Set/CMA1/1_CMA1/' + j + '/DispsInit/' + k,os.getcwd()+'/'+k)
        os.chdir('..')
        os.chdir('..')
        os.chdir('..')
