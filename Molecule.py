# Class for collecting and processing CMA benchmarking data 
# As well as generating LaTeX-formatted SI on a molecule by molecule basis.

import os
import re

import numpy as np
import pandas as pd
pd.set_option("display.max_columns", 15)

class Molecule(object):

    def __init__(self,job):
        info = re.search(r"/(\d*)_.*/(\d*)_(.*)/", job)
        self.ID = f"{info.group(1)}.{info.group(2)}"
        self.name = info.group(3)
        self.geoms = {}
        self.freqs = {}
        self.ted = {}
        self.proj = None
    
    def run(self):
        # ID info
        print('Molecule Class!')
        print("="*25)
        print(f"{self.name + ' ' + self.ID:^25}")
        print("="*25)

        # Geometries
        print("Geometries:")
        print("-----------")
        for geom in self.geoms:
            print(geom)
            for line in self.geoms[geom]:
                print(f"{line[0]:<2}    {float(line[1]):>13.10f}    {float(line[2]):>13.10f}    {float(line[3]):>13.10f}")
            print()

        # Frequencies
        print("Frequencies:")
        print("------------")
        # Creat multiline index
        t = [(head.split()[0], head.split()[1]) for head in self.freqs.keys()] 
        ind = pd.MultiIndex.from_tuples(t)
        df = pd.DataFrame(data=self.freqs)
        df.columns = ind
        print(df.to_string(index=False, float_format="%.2f"))
        print()

        # Nattys
        print("Nattys:")
        print("-------")
        for i, nic in enumerate(self.nics):
            if nic[0] > 1:
                print(f"{i+1:2}    1/\u221a{nic[0]}[", end='')
                for j, term in enumerate(nic[1:]):
                    if j == 0:
                        if term[0] == 1.0:
                            print(term[1], end='')
                        else:
                            print(f"{int(term[0])}{term[1]}", end='')
                        continue
                    if term[0] > 0:
                        sign = "+"
                    else:
                        sign = "-"
                    if abs(term[0]) == 1.0:
                        print(f" {sign} {term[1]}", end='')
                    else:
                        print(f" {sign} {int(abs(term[0]))}{term[1]}", end='')
                print("]")
            else:
                print(f"{i+1:2}    {nic[1][1]}")
        print()
        
        # TED
        # print("TEDs:")
        # print("-----")
        # for ted in self.ted:
        #     print(ted)
        #     print("-"*(13+10*len(self.ted[ted])))
        #     print("     Freq    ", end="")
        #     for freq in self.freqs[f"Natty ({ted})"]:
        #         print(f"{freq:>6.1f}    ", end="")
        #     print("\n"+"-"*(13+10*len(self.ted[ted])))
        #     for i, line in enumerate(self.ted[ted]):
        #         print(f"Mode #{i+1:>3}    ", end="")
        #         for val in line:
        #             print(f"{val:>6.1f}    ", end="")
        #         print()
        #     print("-"*(13+10*len(self.ted[ted]))+"\n")

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
            try: 
                with open(filename, "r") as f:
                    txt = f.readlines()
                self.direc_complete = True 
            except:
                print(f'There is no file named {filename} in this directory')
                self.direc_complete = False
                break
                #self.geoms = None 
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

    def get_nattys(self,combo):
        # Get zmat information
        with open(f"./{combo[0]}/zmat", "r") as f:
            txt = f.readlines()

        zmat = []
        read = False
        for line in txt:
            if "ZMAT begin" in line:
                read = True
                continue
            elif "ZMAT end" in line:
                read = False
                break
            if re.search(r"\d+", line):
                zmat.append(line.split())

        # Format zmat to be more readable
        for i, coord in enumerate(zmat):
            # Bond distances
            if len(coord) == 2:
                zmat[i] = f"r({ coord[0] },{ coord[1] })"
            # Bond angles
            if len(coord) == 3:
                zmat[i] = f"\u03b8({ coord[0] },{ coord[1] },{ coord[2] })"
            if len(coord) == 5:
                # Torsions
                if coord[-1] == "T":
                    zmat[i] = f"\u03c4({ coord[0] },{ coord[1] },{ coord[2] },{ coord[3] })"
                # Out of plane angles
                if coord[-1] == "O":
                    zmat[i] = f"\u03b3({ coord[0] },{ coord[1] },{ coord[2] },{ coord[3] })"
                # Linear bends
                if coord[-1] == "L":
                    zmat[i] = f"\u03b8({ coord[0] },{ coord[1] },{ coord[2] },{ coord[3] })"
                # Linx bends
                if coord[-1] == "Lx":
                    zmat[i] = f"\u03b8x({ coord[0] },{ coord[1] },{ coord[2] },{ coord[3] })"
                # Liny bends
                if coord[-1] == "Ly":
                    zmat[i] = f"\u03b8y({ coord[0] },{ coord[1] },{ coord[2] },{ coord[3] })"
        
        # Combine proj and zmat 
        nics = []
        for coord in self.proj.T:
            tmp = coord/np.abs(coord)[coord>0].min()
            tmp = np.rint(tmp)
            nic = [int(np.sum(np.square(tmp)))]
            for i, coeff in enumerate(tmp):
                if abs(coeff) >= 1:
                    nic.append((coeff, zmat[i]))
            nics.append(nic)

        self.nics = nics


        # Generate NIC representaions
        # for i, coord in enumerate(zmat):
        #     # Bond distances
        #     if len(coord) == 2:
        #         zmat[i] = f"r_{{{ coord[0] },{ coord[1] }}}"
        #     # Bond angles
        #     if len(coord) == 3:
        #         zmat[i] = f"\\theta_{{{ coord[0] },{ coord[1] },{ coord[2] }}}"
        #     if len(coord) == 5:
        #         # Torsions
        #         if coord[-1] == "T":
        #             zmat[i] = f"\\tau_{{{ coord[0] },{ coord[1] },{ coord[2] },{ coord[3] }}}"
        #         # Out of plane angles
        #         if coord[-1] == "O":
        #             zmat[i] = f"\gamma_{{{ coord[0] },{ coord[1] },{ coord[2] },{ coord[3] }}}"
        #         # Linear bends
        #         if coord[-1] == "L":
        #             zmat[i] = f"\\theta_{{{ coord[0] },{ coord[1] },{ coord[2] },{ coord[3] }}}"
        #         # Linx bends
        #         if coord[-1] == "Lx":
        #             zmat[i] = f"\\theta_{{{ coord[0] },{ coord[1] },{ coord[2] },{ coord[3] }}}"
        #         # Liny bends
        #         if coord[-1] == "Ly":
        #             zmat[i] = f"\\theta_{{{ coord[0] },{ coord[1] },{ coord[2] },{ coord[3] }}}"
        # print(*zmat)

