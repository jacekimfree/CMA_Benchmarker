import shutil as sh
import os
import sys

here = os.getcwd()
dirs = os.listdir(here)

for i in dirs:
    if os.path.isdir(os.getcwd() + "/" + i):
        os.chdir(i)
        os.chdir("CCSD_T_TZ")
        print(i)
        # sh.copy(os.getcwd() + '/' + 'zmat',os.getcwd() + '/' + 'zmat_cma1')
        sh.copy(os.getcwd() + '/' + 'zmat_cma1',os.getcwd() + '/' + 'zmat_cma1_Final')
        os.chdir('..')
        os.chdir('..')
