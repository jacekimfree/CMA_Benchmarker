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

pd.set_option("display.max_columns", 15)

# =======================
# Database Specifications
# =======================

# High and low levels of theory
# Available: "CCSD_T_TZ", "CCSD_T_DZ", "B3LYP_6-31G_2df,p_"
h_theory = ["CCSD_T_TZ"]
l_theory = ["CCSD_T_DZ", "B3LYP_6-31G_2df,p_"]
combos = list(product(h_theory,l_theory))

# Coordinates types to use
# Available: "Nattys", "Redundant", "ZMAT" (not yet tho)
coord_type = ["Nattys", "Redundant"]

# Specify paths to grab data from
# Options: '/1_Closed_Shell', '/1_Linear', '/1*', '/2_Open_Shell', '/2*'
paths = ['/2*']

# Various output control statements
n = 0                   # Number of CMA2 corrections (n = 0 -> CMA0)
cma1 = False            # Run CMA1 instead of CMA0
csv = False             # Generate database .csv file
SI = True               # Generate LaTeX SI file
compute_all = False     # Not implemented yet, run calculations for all

if n > 0 and cma1 == True:
    raise RuntimeError("Don't do both CMA1 and CMA2 at the same time yet")

# =====================
# Some useful functions
# =====================

def freq_diff(F1, F2):
    return F1 - F2

def averages(F_diff):
    return np.average(F_diff), np.average(abs(F_diff)) 

def stdev(F_diff):
    return np.std(F_diff)

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
      CCCCCCCC      MMMMMM     MMMMMM           AAAAA
    CCCCCCCCCCCC    MMMMMMM   MMMMMMM          AAAAAAA
   CCCC      CCCC   MMMMMMMM MMMMMMMM         AAAA AAAA
  CCCC              MMMM MMMMMMM MMMM        AAAA   AAAA
  CCCC              MMMM  MMMM   MMMM       AAAA     AAAA
  CCCC              MMMM         MMMM      AAAA       AAAA
  CCCC              MMMM         MMMM     AAAAAAAAAAAAAAAAA 
   CCCC      CCCC   MMMM         MMMM    AAAAAAAAAAAAAAAAAAA
    CCCCCCCCCCCC    MMMM         MMMM   AAAA             AAAA
      CCCCCCCC      MMMM         MMMM  AAAA               AAAA
"""
)

print("Authors: The Dorktor, Nathaniel Kitzpapi, the other one")
print()
print("Combinations of levels of theory (high, low): ", end="")
print(*combos, sep=", ")
print("Coordinate types: ", end="")
print(*coord_type, sep=", ")
print(f"Number of off-diagonal corrections: {n}")
print()

# Grab jobs from each path and order them numerically by ID
for path in paths:
    path_ind.append(len(jobb_list))
    tmp_list = glob.glob(hq + path + "/[1-9]*_*/")
    ind = np.argsort(np.array([int(re.search(r"/\d_.*/(\d*)_.*", name).group(1)) for name in tmp_list]))
    tmp_list = [tmp_list[i] for i in ind]
    jobb_list += tmp_list

# Gives printout for each molecule in jobb_list
print(f"Generating database entries for {len(jobb_list)} jobs:",end="")
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

            # Initialize objects
            mol = Molecule(job)
            basename = return_path_base(job) 
            d = {'Molecule' : None}     # Ensures molecule is first column
            if n > 0: 
                d2 = {'Molecule' : None}

            if i in path_ind:
                print(f"Currently running jobs in {paths[path_ind.index(i)]}\n")
            
            print("////////////////////////////////////////////")
            print(f"//{basename:^40}//")
            print("////////////////////////////////////////////")

            # Run CMA for each combination of theory
            for combo in combos:

                # Grab geometry information 
                mol.get_geoms(combo)
            
                if cma1 == False:
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
                        print(f"  Job                     {mol.name} ({mol.ID})") 
                        print(f"  High level of theory    {combo[0]}")
                        print(f"  Low level of theory     {combo[1]}")
                        print(f"  Coordinate type         {coord}")
                        print("="*50)
                        print()
                    
                        # Skip redundants if linear molecule
                        if 'Linear' in job and coord == 'Redundant':
                            print("What kind of idiot would use redundants for a linear molecule")
                            continue

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
                        execMerger.run(execMerger.options, Proj)
                    
                        # Collect data
                        if coord_type.index(coord) == 0:
                            d[f'Ref ({combo[0]})'] = execMerger.reference_freq
                            
                            # Number the modes
                            d['Molecule'] = [f"{mol.name} ({mol.ID}) mode {i+1}" for i in range(len(execMerger.reference_freq))]
                    
                        if coord == "Nattys":
                            d[f'Natty ({combo[1]})'] = execMerger.Freq_custom
                            d[f'Ref - Nat ({combo[1]})'] = freq_diff(execMerger.reference_freq, execMerger.Freq_custom)
                        if coord == "Redundant":
                            d[f'Red ({combo[1]})'] = execMerger.Freq_redundant
                            d[f'Ref - Red ({combo[1]})'] = freq_diff(execMerger.reference_freq, execMerger.Freq_redundant)
                    
                        # Collect data for CMA2
                        if n > 0:
                            d2[f"Ref {combo[0]}"] = execMerger.reference_freq
                            cma2_data = [] 
                            cma2_dict = execMerger.cma2_dict
                            key = 'cma2_'  + str(execMerger.options.coords) 
                            cma2_data.append(cma2_dict[key]) 
                            print('cma2_data') 
                            print(cma2_data)
                            print(cma2_data[0])
                    
                            if coord == "Nattys":
                                d2[f'Ref - Nat {combo[1]}'] = freq_diff(execMerger.reference_freq, cma2_dict[key])
                            if coord == "Redundant":
                                d2[f'Ref - Red {combo[1]}'] = freq_diff(execMerger.reference_freq, cma2_dict[key])
                    
                        # delete objects so they are forced to reload
                        del execMerger
                        del Merger
                        if coord == "Nattys":
                            del project_obj
                    
                    # end of coord loop

                    # Difference between Natty and Redundant freqs
                    if 'Nattys' in coord_type and 'Redundant' in coord_type and 'Linear' not in job:
                        d[f'Nat - Red {combo[1]}'] = freq_diff(d[f'Natty ({combo[1]})'], d[f'Red ({combo[1]})'])
                    
                    if n > 0:
                        cma2_1_nat_diff = freq_diff(Freq_ref, cma2_data_nat)         
                        cma2_1_red_diff = freq_diff(Freq_ref, cma2_data_red)         
                    
                #======================
                # Mitchell's playground
                #======================

                elif cma1 == True:
                    # move into directory with higher level geom, fc.dat, and Disp directories
                    os.chdir(f"{job}/{combo[0]}/") 
                    
                    # Add your code here
                    print(f"I am in {os.getcwd()} and I can see {os.listdir()}")
                    
                    # In the future we can use this coord loop if we want to 
                    # But for now we will stick with redundants
                    for coord in coord_type:
                        if coord == "Nattys":
                            print("No Nattys :(")
                    
                        else:
                            # Import Merger
                            from Merger import Merger
                            execMerger = Merger()
                            execMerger.options.man_proj = False
                            execMerger.options.coords = coord
                            Proj = None
                            
                            # Run the thing
                            #execMerger.run(execMerger.options, Proj)

                            # Collect the data in dictionary d to add it to the database
                            # e.g. d[f"Ref {combo[0]}"] = execMerger.reference_freq
                            # d[f"CMA1 {combo[1]}"] = execMerger.Freq_redundant

                            del execMerger
                            del Merger

            # end of combo loop

            # Print molecule information
            mol.run()
            
            # Clean up job directory
            if cma1 == False:
                os.remove("zmat")
                os.remove("zmat2")
                os.remove("fc.dat")
                os.remove("fc2.dat")
            sys.path.remove(job)
            del mol

            # Skip dataframe if linear molecule
            if 'Linear' in job:
                continue

            # Add to pandas dataframe
            if csv == True:
                df = pd.DataFrame(data=d)
                # print(df)
                # print()
                frame.append(df)
            
            if n > 0:
                d2 ={'Ref_Redundant' : F2diff,
                     'CMA2_1_el' : cma2_1_red_diff}        
                df2 = pd.DataFrame(data=d2)
                df2.style.apply(highlight_greaterthan, axis =1)
                frame2.append(df2)

        # end of job loop

    # end of def
 
execute()
os.chdir(hq)
if csv == True:
    megaframe = pd.concat(frame)
    megaframe.to_csv('CoordDep.csv', index=False, float_format="%.2f")
    if n > 0:
        megaframe2 = pd.concat(frame2)
        megaframe2.to_csv('CMA2_Convergent.csv', index=False, float_format="%.2f")
#megaframe2.to_excel('CMA2_Convergent.xlsx', index=False)
#writer = pd.ExcelWriter("CMA2_Convergent.xlsx",engine = 'xlsxwriter', encoding = 'utf-8')
#megaframe2.to_excel(writer, sheet_name = 'Sheet1', index=False)
#
#workbook = writer.book
#worksheet = writer.sheets['Sheet1']
#
#workbook.close()



