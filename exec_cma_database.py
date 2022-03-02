import os
import sys
import shutil
import importlib.util
from os import path
import glob
import pandas as pd
import numpy as np
from os.path import dirname, realpath
from runpy import run_path
import copy
import re
from itertools import product

from Molecule import Molecule

# High and low levels of theory, will use all combos
# Available: "CCSD_T_TZ", "CCSD_T_DZ", "B3LYP_6-31G_2df,p_"
h_theory = ["CCSD_T_TZ"]
l_theory = ["CCSD_T_DZ", "B3LYP_6-31G_2df,p_"]
combos = list(product(h_theory, l_theory))

# Coordinates types to use
# Available: "Nattys", "Redundant", "ZMAT" (not yet tho)
coord_type = ["Nattys", "Redundant"]

# Specify paths to grab data from
paths = ['/2_Open_Shell']
#paths = ['/1_Closed_Shell','/2_Open_Shell']

# Number of CMA2 corrections (n=0 -> CMA0)
n = 0

def return_path_base(jobpath):
    name = os.path.basename(os.path.normpath(jobpath)) 
    return name

hq = os.getcwd()
jobb_list = []
path_ind = []

print("CMA Database Generation\n")
print(
"""
      CCCCCCCC        MMMMMM     MMMMMM           AAAAA
    CCCCCCCCCCCC      MMMMMMM   MMMMMMM          AAAAAAA
   CCCC      CCCC     MMMMMMMM MMMMMMMM         AAAA AAAA
  CCCC                MMMM MMMMMMM MMMM        AAAA   AAAA
  CCCC                MMMM  MMMM   MMMM       AAAA     AAAA
  CCCC                MMMM         MMMM      AAAA       AAAA
  CCCC                MMMM         MMMM     AAAAAAAAAAAAAAAAA 
   CCCC      CCCC     MMMM         MMMM    AAAAAAAAAAAAAAAAAAA
    CCCCCCCCCCCC      MMMM         MMMM   AAAA             AAAA
      CCCCCCCC        MMMM         MMMM  AAAA               AAAA
"""
)
print()
print("Combinations of theory to be run (high, low): ", end="")
for combo in combos:
    print("("+combo[0]+", "+combo[1]+")  ", end="")
print("\nCoordinate types: ", end="")
for coord in coord_type:
    print(coord, end="  ")
print()

# Grab jobs from each path and order them numerically by ID
for path in paths:
    tmp_list = glob.glob(hq + path + "/[1-9]*_*/")
    ind = np.argsort(np.array([int(re.search(hq + path + r"/(\d*)_", name).group(1)) for name in tmp_list]))
    tmp_list = [tmp_list[i] for i in ind]
    path_ind.append(len(jobb_list))
    jobb_list += tmp_list

# Gives printout for each molecule in jobb_list
print("Generating database entries for:",end="")
for i, job in enumerate(jobb_list):
    if i in path_ind:
        print("\n\nDirectory: {}".format(paths[path_ind.index(i)]),end="")
        count = 0
    if count % 5 == 0:
        print("\n{0:24}".format(return_path_base(job)),end="")
        count += 1
    else:
        print("{0:24}".format(return_path_base(job)),end="")
        count += 1
print("\n")

compute_all = False

#frame = []
def execute():
    if compute_all:
        print('this feature is not supported in the beta version.')
    else:
        for i, job in enumerate(jobb_list):
            os.chdir(job)
            sys.path.insert(0,job)
            options = None 

            basename = return_path_base(job) 
            print("////////////////////////////////////////////")
            print("//{:^40s}//".format(basename))
            print("////////////////////////////////////////////")
 
            # Run CMA for each combination of theory
            for combo in combos:

                # Copy the necessary files with correct names
                shutil.copyfile(job + combo[1] + "/zmat", job + "zmat")
                shutil.copyfile(job + combo[1] + "/fc.dat", job + "fc.dat")
                shutil.copyfile(job + combo[0] + "/zmat", job + "zmat2")
                shutil.copyfile(job + combo[0] + "/fc.dat", job + "fc2.dat")       

                # Run for each coord type
                for coord in coord_type:

                    # Kept in for its historical significance
                    if coord == "Nattys":
                        print('ReeeeEEEEEeEEEEEEEEEEEEeEeeeeeeeeeeeeeeeeeEEEE')
                    elif coord == "Redundant":
                        print("Catalina wine mixer " + str(i))

                    print()
                    print("="*50)
                    print(" "*17+"Current Analysis")
                    print("-"*50)
                    print("  High level of theory    " + combo[0])
                    print("  Low level of theory     " + combo[1])
                    print("  Coordinate type         " + coord)
                    print("="*50)
                    print()
                
                    from Merger import Merger
                    execMerger = Merger()
                     
                    #Specify options for Merger
                    if coord == "Nattys":
                        execMerger.options.man_proj = True
                        execMerger.options.coords = 'Custom'
            
                        #import the manual_projection module from the specific molecule directory 
                        spec = importlib.util.spec_from_file_location("manual_projection",  job + "/manual_projection.py")
                        foo = importlib.util.module_from_spec(spec)
                        spec.loader.exec_module(foo)
                        project_obj = foo.Projection(None)
                        project_obj.run()
                        Proj = copy.copy(project_obj.Proj)
                    
                    else:
                        execMerger.options.man_proj = False
                        execMerger.options.coords = coord
                        Proj = None

                    execMerger.options.n_cma2 = n

                    # Run CMA
                    output_obj = execMerger.run(execMerger.options, Proj)
                    
                    #Freq_reference = execMerger.reference_freq    
                    #F_nat = execMerger.F_custom 
                    #Freq_nat = execMerger.Freq_custom 
                    
                    #execMolecule = Molecule(natty_obj,redundant_obj,basename)
                    #execMolecule.run()             
                    #F_red = execMerger.F_redundant
                    #Freq_red = execMerger.Freq_redundant 
                    
                    #take differences
                    #F1diff = freq_diff(Freq_reference, Freq_nat)   
                    #F2diff = freq_diff(Freq_reference, Freq_red)   
                    #F1F2diff = freq_diff(Freq_nat, Freq_red)   
                    
                    #dir_of_file = dirname(os.getcwd())
                    #basename = return_path_base(dir_of_file) 
                    #d = {'Molecule': None, 'Ref. Freq': Freq_reference, 'Natty Freq': Freq_nat,'Ref - Nat':F1diff, 'Redundant Freq' : Freq_red, 'Ref - Redundant' :F2diff, 'Nat - Redundant': F1F2diff}        
                    #frame = build_dataframe(d,basename)
                     
                    #print('printing all variables')
                    #print(dir())
                    #delete Merger, execMerger so it is forced to reload all CMA stuff
        
                    del execMerger
                    del Merger
                    if coord == "Nattys":
                        del project_obj
                    
                # end of coord loop
            
            # end of combo loop

            os.remove("zmat")
            os.remove("zmat2")
            os.remove("fc.dat")
            os.remove("fc2.dat")

            sys.path.remove(job)
            # print('printing all variables again')
            # print(dir())
            # os.chdir('../')

        # end of job loop

 
execute()
os.chdir(hq)
#megaframe = pd.concat(frame)
#megaframe.to_csv('chonk.csv', index=False)


