import os
from os import path
import importlib.util

"""
The purpose of this script is open-ended but is intended to do the following:
	1. Generate a database directory tree
	2. Execute CMA scripts and reap data
	3. Upload these results to a Jupyter notebook

"""

#directory_dict = {

#top_layer : 'num_molecules',


l_one = ['1_Closed_Shell', '2_Open_Shell']
l_two = [104, 18]
l_three = ['1_Natural_Internal_Coordinates', '2_Redundant_Internal_Coordinates']
path = os.getcwd()
#def generate():
#    for i, x in enumerate(l_one):
#        for j in range(0, l_two[i]):
#            for k, y in enumerate(l_three):
#                yee = os.path.join(path, x, str(j+1), y)
#                print(yee)
#                #os.sys('mkdir -p' +  '/' + x + '/' + str(j+1) + '/' + y)
#                print(x, j +1, y)
#                os.makedirs(yee)
#    return None

#generate()


l_one_custom = []
l_two_custom = []
l_three_custom = []
print(path)

#yeet = os.path.join(path,l_one[0], '1', l_three[0]) 
#print(yeet)
#print(type(l_one[0]))
#print(type(l_three[0]))

job_list = [
os.path.join(path, l_one[0], '1', l_three[0]),
os.path.join(path, l_one[0], '2', l_three[0]),
os.path.join(path, l_one[1], '2', l_three[0]),
]

print(job_list)
compute_all = False

def execute():
    if compute_all:
        for i, x in enumerate(l_one):
            for j in range(0, l_two[i]):
                for k, y in enumerate(l_three):
                    yee = os.path.join(path, x, str(j+1), y)
                    os.chdir(yee)
                    os.sys('nohup python -u CoordDep.py &')
    else:
        for job in job_list:
            os.chdir(job)
            spec = importlib.util.spec_from_file_location("module.name", "/path/to/file.py")
            foo = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(foo)
            foo.MyClass()            
            suc.some_func()
    return None
execute()     
                   
