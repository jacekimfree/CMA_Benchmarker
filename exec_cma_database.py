import os
import sys
import importlib.util
from os import path
import glob
import pandas as pd
import numpy as np
from os.path import dirname, realpath
from runpy import run_path
import copy
import re
from Molecule import Molecule

#job_prefixes = [18] #[80] #[21,58,59,57,65,92,93,55,56,11,13,14,3,53,54,60,63,20,62,61,49,50,16,15,2,6,22,24,25,29,30,31,33,67,23,101,35,36,8,96,103,26,27,28] #17,18,19, 47
#jobb_list = []
#for jobs in job_prefixes:
#    prefix = glob.glob(path + "/" + str(jobs) + "_*")[0]
#    j = [os.path.join(prefix, '1_Natural_Internal_Coordinates'), os.path.join(prefix, '2_Redundant_Internal_Coordinates')]
#    jobb_list.append(j)

n = 2

hq = os.getcwd()
path = os.getcwd() + '/1_Closed_Shell'

jobb_list = glob.glob(path + "/*")
jobb_list = filter(lambda x: not re.search("/hide", x), jobb_list)
print(jobb_list)

compute_all = False

frame = []
def return_path_base(jobpath):
    name = os.path.basename(os.path.normpath(jobpath)) 
    return name

def execute():
    if compute_all:
        print('this feature is not supported in the beta version.')
    else:
        for i, job in enumerate(jobb_list):
            job = glob.glob(job + "/1_Natural_Internal_Coordinates")
            os.chdir(job[0])
            print('ReeeeEEEEEeEEEEEEEEEEEEeEeeeeeeeeeeeeeeeeeEEEE')
            print(os.getcwd()) 
            sys.path.insert(0,job[0])
            options = None 
                       
            #import the manual_projection module from the specific molecule directory 

            spec = importlib.util.spec_from_file_location("manual_projection",  job[0] + "/manual_projection.py")
            foo = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(foo)
            project_obj = foo.Projection(None)
            
            dir_of_file = dirname(os.getcwd())
            basename = return_path_base(dir_of_file) 
            print("////////////////////////////////////////////")
            print("//{:^40s}//".format(basename))
            print("////////////////////////////////////////////")
 
            #import manual_projection as manual_projection 
            project_obj.run()
            
            from Merger import Merger
            execMerger = Merger()
             
            #Specify options for Nat_int_coords

            execMerger.options.man_proj = True
            execMerger.options.coords = 'Custom'
            execMerger.options.n_cma2 = n 
            natty_obj = execMerger.run(execMerger.options, copy.copy(project_obj.Proj))
            
            #Freq_reference = execMerger.reference_freq    
            #F_nat = execMerger.F_custom 
            #Freq_nat = execMerger.Freq_custom 
            
            #move to the redundant internal coordinates directory
            os.chdir('../') 
            print('directory')
            print(os.getcwd())
            #os.chdir(job[1])     
            os.chdir('2_Redundant_Internal_Coordinates')
            sys.path.remove(job[0])
            sys.path.insert(0, job[0])
            
            #Redo everything, but with redundant internal coordinates
            
            print("Catalina wine mixer " + str(i))
            from Merger import Merger
            execMerger = Merger()
            #Specify options for redunant_int_coords

            execMerger.options.man_proj = False
            execMerger.options.coords = 'Redundant'
            execMerger.options.n_cma2 = n 
            redundant_obj = execMerger.run(execMerger.options, None)
            execMolecule = Molecule(natty_obj,redundant_obj,basename)
            execMolecule.run()             
            #F_red = execMerger.F_redundant
            #Freq_red = execMerger.Freq_redundant 
            
            #take differences
            #F1diff = freq_diff(Freq_reference, Freq_nat)   
            #F2diff = freq_diff(Freq_reference, Freq_red)   
            #F1F2diff = freq_diff(Freq_nat, Freq_red)   
            
            dir_of_file = dirname(os.getcwd())
            basename = return_path_base(dir_of_file) 
            #d = {'Molecule': None, 'Ref. Freq': Freq_reference, 'Natty Freq': Freq_nat,'Ref - Nat':F1diff, 'Redundant Freq' : Freq_red, 'Ref - Redundant' :F2diff, 'Nat - Redundant': F1F2diff}        
            #frame = build_dataframe(d,basename)
             
            print('printing all variables')
            print(dir())
            #delete Merger, execMerger so it is forced to reload all CMA stuff
            sys.path.remove(job[0])

            del execMerger
            del Merger
            del project_obj
            print('printing all variables again')
            print(dir())
            os.chdir('../')

 

execute()
os.chdir(hq)
#megaframe = pd.concat(frame)
#megaframe.to_csv('chonk.csv', index=False)


