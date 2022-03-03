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

# =======================
# Database Specifications
# =======================

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

compute_all = False

# =====================
# Some useful functions
# =====================

def freq_diff(F1, F2):
    return F1 - F2

def averages(F_diff):
    return np.average(F_diff), np.average(abs(F_diff)) 
def stdev(F_diff):
    return np.std(F_diff)

def return_path_base(jobpath):
    name = os.path.basename(os.path.normpath(jobpath)) 
    return name

def highlight_greaterthan(x):
    print('this is x.Ref__Redundant')
    print(x.Ref_Redundant) 
    #return np.where((x.F2diff > 0.5), props, '')
    if abs(x.Ref_Redundant) > 0.5:
        print('oi!!!!') 
        return ['background-color:yellow']
    else:
         return ['background-color:white']

def build_dataframe(d,basename):
    df = pd.DataFrame(data=d)
    df_i = df.to_csv('out.csv',index=False)
    # adds a space in between each molecule. 
    df.loc[df.shape[0]] = [None, None, None, None, None, None, None] 
    df.loc[df.shape[0]] = [None, None, None, None, None, None, None] 
    df.loc[df.shape[0]] = [None, None, None, None, None, None, None] 

    #df.loc[df.shape[0]] = [None, None, 'Signed Avg.', averages(F1diff)[0], 'Signed Avg.', averages(F1diff)[0], averages(F1F2diff)[0]]
    #df.loc[df.shape[0]] = [None, None, 'Average   .', averages(F1diff)[1], 'Average    ', averages(F1diff)[1], averages(F1F2diff)[1]]
    #df.loc[df.shape[0]] = [None, None, 'Std. Dev  .', stdev(F1diff),    'Std. Dev.  ', stdev(F1diff), stdev(F1F2diff)]

    df.loc[0,('Molecule')] = str(basename) 
    #df_i: dataframe, individual for each child directory
    frame.append(df)
    return frame, df_i 

def return_path_base(jobpath):
    name = os.path.basename(os.path.normpath(jobpath)) 
    return name


# ====================================
# Print out information about database
# ====================================

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
        print("\n{0:20}".format(return_path_base(job)),end="")
        count += 1
    else:
        print("{0:20}".format(return_path_base(job)),end="")
        count += 1
print("\n")

# Initialize frames for pandas database
frame = []
if n > 0:
    frame2 = []


# ============
# Do the thing
# ============

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

                # Lets us know we need to save reference data for current theory combo
                Freq_ref = None   

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

                    # Collect data
                    if Freq_ref == None:
                        Freq_ref = execMerger.reference_freq    
                    if coord == "Nattys":
                        Freq_nat = execMerger.Freq_custom
                    if coord == "Redundant":
                        Freq_red = execMerger.Freq_redundant

                    # Collect data for CMA2
                    if n > 0:
                        cma2_data = [] 
                        cma2_dict = execMerger.cma2_dict
                        for z in range(0,n):
                            key = 'cma2_'  + str(execMerger.options.coords) 
                            cma2_data.append(cma2_dict[key]) 
                        print('cma2_data') 
                        print(cma2_data)
                        print(cma2_data[0])

                        if coord == "Nattys":
                            cma2_data_nat = cma2_data[0]
                        if coord == "Redundant":
                            cma2_data_red = cma2_data[0]

                    # delete objects so they are forced to reload
                    del execMerger
                    del Merger
                    if coord == "Nattys":
                        del project_obj
                    
                # end of coord loop

                # take differences
                F1diff = freq_diff(Freq_ref, Freq_nat)   
                F2diff = freq_diff(Freq_ref, Freq_red)   
                F1F2diff = freq_diff(Freq_nat, Freq_red)   
                if n > 0:
                    cma2_1_nat_diff = freq_diff(Freq_ref, cma2_data_nat)         
                    cma2_1_red_diff = freq_diff(Freq_ref, cma2_data_red)         
            
                dir_of_file = dirname(os.getcwd())
                basename = return_path_base(dir_of_file) 
                d = {'Molecule': basename, 'Ref. Freq': Freq_ref, 'Natty Freq': Freq_nat,
                     'Ref - Nat':F1diff, 'Redundant Freq' : Freq_red, 'Ref - Redundant' :F2diff, 'Nat - Redundant': F1F2diff}
                df = pd.DataFrame(data=d)
                frame.append(df)
                
                if n > 0:
                    d2 ={'Ref_Redundant' :F2diff,'CMA2_1_el' : cma2_1_red_diff}        
                    df2 = pd.DataFrame(data=d2)
                    df2.style.apply(highlight_greaterthan, axis =1)
                    frame2.append(df2)

            # end of combo loop

            # Clean up job directory
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
megaframe = pd.concat(frame)
megaframe.to_csv('CoordDep.csv', index=False)
if n > 0:
    megaframe2 = pd.concat(frame2)
    megaframe2.to_csv('CMA2_Convergent.csv', index=False)
#megaframe2.to_excel('CMA2_Convergent.xlsx', index=False)
#writer = pd.ExcelWriter("CMA2_Convergent.xlsx",engine = 'xlsxwriter', encoding = 'utf-8')
#megaframe2.to_excel(writer, sheet_name = 'Sheet1', index=False)
#
#workbook = writer.book
#worksheet = writer.sheets['Sheet1']
#
#workbook.close()



