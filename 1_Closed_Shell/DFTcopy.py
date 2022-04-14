# Quick copy script for grabing the DFT data from original files in hm97506
import os
import shutil as sh
import re
from glob import glob 

dirs = glob("*/")

for d in dirs:
    num = re.search("(\d*)_.*",d).group(1)
    path = f"/home/vulcan/hm97506/CMA/G2_benchmark/b3lyp/1_ClosedShell/{num}/2_CMA/"

    with open(path+"zmat","r") as f:
        oldzmat = f.readlines()
    
    read = False
    dft_cart = []
    for line in oldzmat:
        if "cart begin" in line:
            read = True
            continue
        elif "---" in line:
            read = False
            break
        if read:
            dft_cart.append(line)

    with open(d+"CCSD_T_TZ/zmat","r") as f:
        tz = f.readlines()

    zmat = []
    for line in tz:
        if "ZMAT begin" in line:
            read = True
            continue
        elif "ZMAT end" in line:
            read = False
            break
        if read:
            zmat.append(line)

    with open(d+"B3LYP_6-31G_2df,p_/zmat","w") as w:
        w.write("ZMAT begin\n")
        w.writelines(zmat)
        w.write("ZMAT end\n")
        w.write("\ncart begin\n")
        w.writelines(dft_cart)
        w.write("cart end")

    sh.copyfile(path+"fc.dat",d+"B3LYP_6-31G_2df,p_/fc.dat")
