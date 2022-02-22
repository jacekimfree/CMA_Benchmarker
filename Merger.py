import os
import re
import sys
import shutil
import subprocess
import time
import numpy as np
from numpy import linalg as LA
from numpy.linalg import inv
from scipy.linalg import fractional_matrix_power
from scipy.linalg import block_diag
from concordantmodes.algorithm import Algorithm
from concordantmodes.directory_tree import DirectoryTree
from concordantmodes.f_convert import FcConv
from concordantmodes.f_read import FcRead
from concordantmodes.force_constant import ForceConstant
from concordantmodes.gf_method import GFMethod
from concordantmodes.g_matrix import GMatrix
from concordantmodes.int2cart import Int2Cart
from concordantmodes.options import Options
from concordantmodes.reap import Reap
from concordantmodes.s_vectors import SVectors
from concordantmodes.submit import Submit
from concordantmodes.ted import TED
from concordantmodes.trans_disp import TransDisp
from concordantmodes.vulcan_template import VulcanTemplate
from concordantmodes.zmat import Zmat
import copy
from fractions import Fraction
#import manual_projection
#manual_projection

class Merger(object):

    def __init__(self):
        #print("nothing to init")
        options_kwargs = {
            "queue": "gen4.q,gen6.q,gen5.q",
            "program": "molpro@2010.1.67+mpi",
            "energy_regex": r"\(T\) energy\s+(\-\d+\.\d+)",
            'energy_regex' : r"\s*\!CCSD\(T\) total energy\s+(-\d+\.\d+)",
            "cart_insert": 9,
            "calc" : False,
            "success_regex": r"Variable memory released",
        }
        options_obj = Options(**options_kwargs)
        self.options = options_obj 
    #function that returns diagonal fc matrix + n-largest off-diagonal elements
    def run(self,opts, Proj):

        print("You have imported the merger script!")
        
        #options_kwargs = {
        #    "queue": "gen4.q,gen6.q,gen5.q",
        #    "program": "molpro@2010.1.67+mpi",
        #    "energy_regex": r"\(T\) energy\s+(\-\d+\.\d+)",
        #    'energy_regex' : r"\s*\!CCSD\(T\) total energy\s+(-\d+\.\d+)",
        #    "cart_insert": 9,
        #    #"man_proj" : True,
        #    "calc" : False,
        #    #"coords": "Custom",
        #    "success_regex": r"Variable memory released",
        #}
        #options_obj = Options(**options_kwargs)
        self.Proj = Proj 
        
        options = opts 
        #options = options_obj
        
        rootdir = os.getcwd()
        zmat_obj = Zmat(options)
        zmat_obj.run()
                
        # Build manual projection matrix here.
        np.set_printoptions(edgeitems=60,linewidth=1000)
        #if options.man_proj: 
        #    project_obj = manual_projection.Projection(options)
        #    project_obj.run()
        #    
        #    
        #    self.Proj = project_obj.Proj
        #    print('why is this not reloading')
        #    print(self.Proj)
        #else:
        #    self.Proj = None 
        # Compute the initial s-vectors
        #frac_rep = copy.copy(Proj)
        ##np.set_printoptions(formatter={'all':lambda x: str(fractions.Fraction(x).limit_denominator())})
        ##print(A_inv)
        #yeet = np.vectorize(Fraction)(frac_rep)
        #print('trying to print a dang fraction')
        #print(yeet)
        s_vec = SVectors(
            zmat_obj, options, zmat_obj.variable_dictionary_init
        )
        print('this is proj, check for this when redundants executed', self.Proj)
        s_vec.run(zmat_obj.cartesians_init, True, proj=self.Proj)
                
        TED_obj = TED(s_vec.proj, zmat_obj)
        
        g_mat = GMatrix(zmat_obj, s_vec, options)
        g_mat.run()
        
        G = g_mat.G.copy()
        Gdz = G.copy()
        init_bool = False
        if os.path.exists(rootdir + "/fc.dat"):
            f_read_obj = FcRead("fc.dat")
        elif os.path.exists(rootdir + "/FCMFINAL"):
            f_read_obj = FcRead("FCMFINAL")
        else:
            raise RuntimeError
        
        if not init_bool:
            f_read_obj.run()
            f_conv_obj = FcConv(
                f_read_obj.fc_mat,
                s_vec,
                zmat_obj,
                "internal",
                False,
                TED_obj,
                options.units,
            )
            f_conv_obj.run()
            F = f_conv_obj.F
        else:
            F = fc_init.FC
        
        if options.coords != "ZMAT" and not init_bool:
            F = np.dot(TED_obj.proj.T, np.dot(F, TED_obj.proj))
        if options.coords != "ZMAT":
            g_mat.G = np.dot(TED_obj.proj.T, np.dot(g_mat.G, TED_obj.proj))
        
        TED_obj.run(np.eye(TED_obj.proj.shape[1]),np.zeros(TED_obj.proj.shape[1]))
                
        print("Initial Frequencies:")
        init_GF = GFMethod(
            g_mat.G.copy(),
            F.copy(),
            options.tol,
            options.proj_tol,
            zmat_obj,
            TED_obj,
            False
        )
        init_GF.run()
        
        # Now for the TED check.
        G = np.dot(np.dot(LA.inv(init_GF.L), g_mat.G), LA.inv(init_GF.L).T)
        G[np.abs(G) < options.tol] = 0
        F = np.dot(np.dot(init_GF.L.T, F), init_GF.L)
        F[np.abs(F) < options.tol] = 0
        
        print("TED Frequencies:")
        TED_GF = GFMethod(
            G,
            F,
            options.tol,
            options.proj_tol,
            zmat_obj,
            TED_obj,
            False
        )
        TED_GF.run()
        
        proj_tol = 1.0e-3
        eig_inv = inv(init_GF.L)  # (Normal modes (Q) x Sym internals (S) )
        for i in range(len(eig_inv)):
            eig_inv[i] = eig_inv[i] / LA.norm(eig_inv[i])
            eig_inv[i][
                np.abs(eig_inv[i]) < np.max(np.abs(eig_inv[i])) * proj_tol
            ] = 0
        
        print("Everything before this statement has been crosschecked with merger/coordep") 
        # Now run the TZ force constant transformation
        
        zmat_obj2 = Zmat(options)
        zmat_obj2.run(zmat_name="zmat2")
        
        #print('final proj type')
        #print(self.Proj)
        #print('now printing ted obj projection thingy')
        #print(TED_obj.proj)
        #print('testing man_proj')
        ### after allowing options.man_proj to = True, I match the redundant int coord script's frequencies exactly. This doesn't make sense
        print(options.man_proj)
        options.man_proj = True
        #s_vec.run(zmat_obj2.cartesians_init, True, proj=Proj)
        
        
        s_vec = SVectors(
            zmat_obj2, options, zmat_obj2.variable_dictionary_init
        )
        s_vec.run(zmat_obj2.cartesians_init, True, proj=TED_obj.proj)
                
        TED_obj = TED(s_vec.proj, zmat_obj2)
                
        g_mat = GMatrix(zmat_obj2, s_vec, options)
        g_mat.run()
        
        G = g_mat.G.copy()
        Gtz = G.copy()
        init_bool = False
        if os.path.exists(rootdir + "/fc2.dat"):
            f_read_obj = FcRead("fc2.dat")
        elif os.path.exists(rootdir + "/FCMFINAL2"):
            f_read_obj = FcRead("FCMFINAL2")
        else:
            raise RuntimeError
        
        if not init_bool:
            f_read_obj.run()
            f_conv_obj = FcConv(
                f_read_obj.fc_mat,
                s_vec,
                zmat_obj2,
                "internal",
                False,
                TED_obj,
                options.units,
            )
            f_conv_obj.run()
            F = f_conv_obj.F
        else:
            F = fc_init.FC
        # redundant basis 
        G = np.dot(np.dot(TED_obj.proj.T,G),TED_obj.proj)
        #G = np.dot(np.dot(eig_inv, G), eig_inv.T)
        # G[np.abs(G) < options.tol] = 0
        F = np.dot(np.dot(TED_obj.proj.T,F),TED_obj.proj)
        #F = np.dot(np.dot(inv(eig_inv).T, F), inv(eig_inv))
        #We shouldn't be zeroing out parts of the FC matrix just because they're small
        #F[np.abs(F) < options.tol] = 0
         
        init_GF = GFMethod(
            G,
            F,
            options.tol,
            options.proj_tol,
            zmat_obj2,
            TED_obj,
            False
        )
        init_GF.run()
        #print('frequencies???')
        self.reference_freq = init_GF.freq 
        #for diagnostic 
        if options.coords == 'Redundants':
            L_B = init_GF.L
        elif options.coords == 'Custom':
            L_A = init_GF.L
         
        #L tensor tempering to normal mode basis
        G = np.dot(np.dot(eig_inv, G), eig_inv.T)
        F = np.dot(np.dot(inv(eig_inv).T, F), inv(eig_inv))
        
        def n_largest(n, FC):
            indexes = []
            upper_triang = np.triu(FC,n)
            print("REEEEEEEEEEEEeeeeeeeeeEEEEEEEEEEE**************")
            for i in range(0,n):
                fc_cma2 = np.where(upper_triang == upper_triang.max())
                index = [fc_cma2[0][0], fc_cma2[1][0]]
                indexes.append(index)
                upper_triang[index[0],index[1]] = 0
            print(indexes)
            return indexes
        
        print("Full Force constant matrix in lower level normal mode basis:")
        print(F)
        if options.coords == 'Redundant':
            self.F_redundant = F
            #self.F_redundant_cma2IDX = n_largest(2, np.abs(copy.copy(self.F_redundant)))
        elif options.coords == 'Custom':
            self.F_custom = F
            #self.F_custom_cma2IDX = n_largest(2, np.abs(copy.copy(self.F_custom)))
        elif options.coords == 'ZMAT' :
            self.F_zmat = F
        else:
            pass
        F = np.diag(np.diag(F))
        
        print("Diagonal Force constant matrix in lower level normal mode basis:")
        print(F)
        #print('G matrix requested by mitchell')
        #print(G) 
        init_GF = GFMethod(
            G,
            F,
            options.tol,
            options.proj_tol,
            zmat_obj2,
            TED_obj,
            False
        )
        
        init_GF.run()
        if options.coords == 'Redundant':
            self.Freq_redundant = init_GF.freq
        elif options.coords == 'Custom':
            self.Freq_custom = init_GF.freq
        elif options.coords == 'ZMAT' :
            self.Freq_zmat = init_GF.freq
        else:
            pass
        print('merger needs to delete stuff as well')
        print(dir())
        #everything below this line pertains to CMA2 off-diag elements being included in the GF matrix computation, its not an optimal setup, obviously
        #plz fix it :( 
        self.cma2_freq = [] 
        if self.options.n_cma2 > 0:
            for index in range(0, self.options.n_cma2):
                if options.coords == 'Redundant':
                    extras = n_largest(2, np.abs(copy.copy(self.F_redundant)))
                    F[index] = self.F_redundant[index]     
                elif options.coords == 'Custom':
                    extras = n_largest(2, np.abs(copy.copy(self.F_custom)))
                    F[index] = self.F_custom[index]     
                elif options.coords == 'ZMAT' :
                    extras = n_largest(2, np.abs(copy.copy(self.F_zmat)))
                    F[index] = self.F_zmat[index]     
                else:
                    pass
                print('Time for some off-diags')
                init_GF = GFMethod(
                    G,
                    F,
                    options.tol,
                    options.proj_tol,
                    zmat_obj2,
                    TED_obj,
                    False
                )
                init_GF.run()
                print('CMA2 including ' + str(index + 1) + ' off-diagonal elements for ' + str(options.coords) + ' coordinates')
                print(init_GF.freq) 

        else:
            print('CMA0 it is') 
