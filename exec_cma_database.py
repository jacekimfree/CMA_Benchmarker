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

#n = 2
#
#hq = os.getcwd()
#path = os.getcwd() + '/1_Closed_Shell'
#
#jobb_list = glob.glob(path + "/*")
#jobb_list = filter(lambda x: not re.search("/hide", x), jobb_list)
#print(jobb_list)

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



#Specify paths to grab data from
#paths = ['/2_Open_Shell']
paths = ['/1_Closed_Shell','/2_Open_Shell']

# Number of CMA2 corrections (n=0 -> CMA0)
n = 2

hq = os.getcwd()
jobb_list = []
path_ind = []

# Grab jobs from each path and order them numerically by ID
for path in paths:
    tmp_list = glob.glob(hq + path + "/[1-9]*_*/")
    ind = np.argsort(np.array([int(re.search(hq + path + r"/(\d*)_", name).group(1)) for name in tmp_list]))
    tmp_list = [tmp_list[i] for i in ind]
    path_ind.append(len(tmp_list) + len(jobb_list))
    jobb_list += tmp_list

# Gives printout for each molecule run
print("Generating database entries for:",end="")
for i, job in enumerate(jobb_list):
    if i == 0 or i in path_ind:
        print("\n\nDirectory: {}".format(paths.pop(0)),end="")
        count = 0
    if count % 5 == 0:
        print("\n{0:20}".format(os.path.basename(job[:-1])),end="")
        count += 1
    else:
        print("{0:20}".format(os.path.basename(job[:-1])),end="")
        count += 1
print("\n")


compute_all = False

def return_path_base(jobpath):
    name = os.path.basename(os.path.normpath(jobpath)) 
    return name

frame = []
frame2 = []
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
            
            Freq_reference = execMerger.reference_freq    
            #F_nat = execMerger.F_custom 
            Freq_nat = execMerger.Freq_custom 
            cma2_data = [] 
            cma2_dict = execMerger.cma2_dict
            for z in range(0,n):
                key = 'cma2_'  + str(execMerger.options.coords) 
                #print(key)
                cma2_data.append(cma2_dict[key]) 
            print('cma2_data') 
            print(cma2_data)
            print(cma2_data[0])
            cma2_1data_nat = cma2_data[0] 
            #cma2_2data_nat = cma2_data[1] 

  
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
            Freq_red = execMerger.Freq_redundant 

            cma2_data = [] 
            cma2_dict = execMerger.cma2_dict
            for z in range(0,n):
                key = 'cma2_'  + str(execMerger.options.coords) 
                #print(key)
                cma2_data.append(cma2_dict[key]) 
            print('cma2_data') 
            cma2_1data_red = cma2_data[0] 
            #cma2_2data_red = cma2_data[1] 


            
            #take differences
            F1diff = freq_diff(Freq_reference, Freq_nat)   
            F2diff = freq_diff(Freq_reference, Freq_red)   
            F1F2diff = freq_diff(Freq_nat, Freq_red)   
            cma2_1_nat_diff = freq_diff(Freq_reference, cma2_1data_nat)         
            #cma2_2_nat_diff = freq_diff(Freq_reference, cma2_2data_nat)         
            cma2_1_red_diff = freq_diff(Freq_reference, cma2_1data_red)         
            #cma2_2_red_diff = freq_diff(Freq_reference, cma2_2data_red)         

 
            dir_of_file = dirname(os.getcwd())
            basename = return_path_base(dir_of_file) 
            d = {'Molecule': basename, 'Ref. Freq': Freq_reference, 'Natty Freq': Freq_nat,'Ref - Nat':F1diff, 'Redundant Freq' : Freq_red, 'Ref - Redundant' :F2diff, 'Nat - Redundant': F1F2diff}
            d2 ={'Ref_Redundant' :F2diff,'CMA2_1_el' : cma2_1_red_diff}        
            #frame = build_dataframe(d,basename)

            df = pd.DataFrame(data=d)
            df2 = pd.DataFrame(data=d2)
            #print(df2.CMA2_1_el)
            df2.style.apply(highlight_greaterthan, axis =1)
            #df2.style.apply(lambda x: ['background:yellow' if abs(x) > 0.5 else 'background:white' for x in df2.F2diff]) 
            #df2.apply(highlight_greaterthan, props = 'color:white;background-color:yellow',axis=1)
            #df_i = df.to_csv('out.csv',index=False)
            # adds a space in between each molecule. 
            #df.loc[df.shape[0]] = [None, None, None, None, None, None, None] 
            #df.loc[df.shape[0]] = [None, None, None, None, None, None, None] 
            #df.loc[df.shape[0]] = [None, None, None, None, None, None, None] 

            #df.loc[df.shape[0]] = [None, None, 'Signed Avg.', averages(F1diff)[0], 'Signed Avg.', averages(F1diff)[0], averages(F1F2diff)[0]]
            #df.loc[df.shape[0]] = [None, None, 'Average   .', averages(F1diff)[1], 'Average    ', averages(F1diff)[1], averages(F1F2diff)[1]]
            #df.loc[df.shape[0]] = [None, None, 'Std. Dev  .', stdev(F1diff),    'Std. Dev.  ', stdev(F1diff), stdev(F1F2diff)]

            #df.loc[0,('Molecule')] = str(basename) 
            
            #df_i: dataframe, individual for each child directory
            frame.append(df)
            frame2.append(df2)

             
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
megaframe = pd.concat(frame)
megaframe2 = pd.concat(frame2)
megaframe.to_csv('CoordDep.csv', index=False)
megaframe2.to_csv('CMA2_Convergent.csv', index=False)
#megaframe2.to_excel('CMA2_Convergent.xlsx', index=False)
#writer = pd.ExcelWriter("CMA2_Convergent.xlsx",engine = 'xlsxwriter', encoding = 'utf-8')
#megaframe2.to_excel(writer, sheet_name = 'Sheet1', index=False)
#
#workbook = writer.book
#worksheet = writer.sheets['Sheet1']
#
#workbook.close()



