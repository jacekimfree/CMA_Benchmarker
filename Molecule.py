# Class for collecting and processing CMA benchmarking data 
# As well as generating LaTeX-formatted SI on a molecule by molecule basis.

import os
import re

import pandas as pd

class Molecule(object):

    def __init__(self,job):
        info = re.search(r"/(\d*)_.*/(\d*)_(.*)/", job)
        self.ID = f"{info.group(1)}.{info.group(2)}"
        self.name = info.group(3)
        self.geoms = {}
        self.freqs = {}
        #self.ted = 
        #self.proj_mat = 
    
    def run(self):
        # ID info
        print('Molecule Class!')
        print("="*25)
        print(f"{self.name + ' ' + self.ID:^25}")
        print("="*25)

        # Geometries
        print("Geometries:")
        for geom in self.geoms:
            print()
            print(geom)
            for line in self.geoms[geom]:
                print(f"{line[0]:<2}    {float(line[1]):>13.10f}    {float(line[2]):>13.10f}    {float(line[3]):>13.10f}")
        print()

        # Frequencies
        print("Frequencies:")
        print()
        print(pd.DataFrame(data=self.freqs).to_string(index=False))
        print()

    def get_geoms(self, combo):
        cma1 = False
        if combo[0] not in os.listdir():
            cma1 = True
        for lvl in combo:
            if lvl in self.geoms.keys():
                if cma1 == True:
                    break
                else:
                    continue
            if cma1 == True:
                filename = f"./zmat"
            else:
                filename = f"./{lvl}/zmat"
            with open(filename, "r") as f:
                txt = f.readlines()
            cart = []
            read = False
            for line in txt:
                if 'cart begin' in line:
                    read = True
                    continue
                elif 'cart end' in line:
                    break

                if read == True:
                    cart.append(line.split())
            self.geoms[lvl] = cart
            if cma1 == True:
                break

        
