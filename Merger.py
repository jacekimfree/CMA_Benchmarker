import os
import re
import sys
import shutil
import json
import subprocess
import time
import numpy as np
import copy

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
from concordantmodes.g_read import GrRead
from concordantmodes.int2cart import Int2Cart
from concordantmodes.options import Options
from concordantmodes.reap import Reap
from concordantmodes.s_vectors import SVectors
from concordantmodes.submit import Submit
from concordantmodes.ted import TED
from concordantmodes.transf_disp import TransfDisp
from concordantmodes.vulcan_template import VulcanTemplate
from concordantmodes.zmat import Zmat

import copy
from fractions import Fraction

class Merger(object):

    def __init__(self, cma1_path=None):
        #print("nothing to init")
        options_kwargs = {
            # "queue": "gen4.q,gen6.q",
            "queue": "gen4.q,gen6.q,gen5.q",
            "program": "molpro@2010.1.67+mpi",
            "energy_regex": r"\(T\) energy\s+(\-\d+\.\d+)",
            'energy_regex' : r"\s*\!CCSD\(T\) total energy\s+(-\d+\.\d+)",
            "cart_insert": 9,
            "calc" : False,
            "calc_init" : False,
            "success_regex": r"Variable memory released",
            # "reduced_disp" : True,
            # "disp" : 1.0
        }
        options_obj = Options(**options_kwargs)
        self.options = options_obj 
        self.cma1_path = cma1_path
    #function that returns diagonal fc matrix + n-largest off-diagonal elements
    def run(self, opts, Proj, energy_regex=None, success_regex=None, cma1_coord=None, sym_sort=None, xi_tol=[], coord_type_init="internal"):
    # def run(self, opts, Proj, energy_regex=None, success_regex=None, cma1_coord=None):
        
        self.coord_type_init = coord_type_init
        
        print("You have imported the merger script!")
        
        self.Proj = Proj
        
        self.lone_wolves = np.array([])
        if len(sym_sort) > 1:
            flat_sym_sort = np.array([])
            for i in range(len(sym_sort)):
                flat_sym_sort = np.append(flat_sym_sort,sym_sort[i])
                if len(sym_sort[i]) == 1:
                    self.lone_wolves = np.append(self.lone_wolves,sym_sort[i])
            flat_sym_sort = flat_sym_sort.astype(int)
            self.lone_wolves = self.lone_wolves.astype(int)
            # self.Proj = self.Proj[:,flat_sym_sort]
        
        if len(self.lone_wolves):
            print("Exclude frequencies with these irreps from the stats: ")
            print(self.lone_wolves)

        options = opts 
        # options.disp = 0.02
        #options = options_obj
        # options.cart_insert_init = 26
        # options.cart_insert_init = 24
        # options.cart_insert_init = 9
        # options.cart_insert_init = 5
        rootdir = os.getcwd()
        zmat_obj = Zmat(options)
        zmat_obj.run()

        np.set_printoptions(edgeitems=60,linewidth=1000)
        
        # Compute the initial s-vectors
        s_vec = SVectors(
            zmat_obj, options, zmat_obj.variable_dictionary_init
        )
        if len(np.shape(self.Proj)) > 2:
            print('this is proj that has been manually sorted by symmetry irrep')
        else:
            print('this is proj, check for this when redundants executed')
            print(self.Proj)
        s_vec.run(zmat_obj.cartesians_init, True, proj=self.Proj, second_order=options.second_order)
                
        TED_obj = TED(s_vec.proj, zmat_obj)
        print("TED PROJ:")
        print(TED_obj.proj)
        
        g_mat = GMatrix(zmat_obj, s_vec, options)
        g_mat.run()
        
        G = g_mat.G.copy()
        
        if os.path.exists(rootdir + "/fc.grad"):
            print('FC GRAD EXISTS')
            # raise RuntimeError
            g_read_obj = GrRead("fc.grad")
            g_read_obj.run(zmat_obj.cartesians_init)
            print(zmat_obj.cartesians_init)

        init_bool = False
        if os.path.exists(rootdir + "/fc.dat"):
            f_read_obj = FcRead("fc.dat")
        elif os.path.exists(rootdir + "/FCMFINAL"):
            f_read_obj = FcRead("FCMFINAL")
        else:
            init_bool = True
            if cma1_coord == None:
                print("You need to specify the cma1_coord variable for this feature. Check execMerger.run()")
                raise RuntimeError

            os.chdir(os.getcwd() + self.cma1_path)
            if os.path.exists(os.getcwd() + "/fc_int_"+cma1_coord+".dat") and not options.second_order:
                if os.path.exists(os.getcwd()+'/DispsInit'):
                    shutil.rmtree("DispsInit")
                f_read_obj = FcRead("fc_int_"+cma1_coord+".dat")
                f_read_obj.run()
                fc_init = ForceConstant(
                    None,
                    [],
                    [],
                    0,
                    options,
                    [],
                )
                fc_init.FC =  f_read_obj.fc_mat
                os.chdir('..')
                os.chdir('..')
            elif os.path.exists(os.getcwd() + "/fc_cart.dat"):
                if os.path.exists(os.getcwd()+'/DispsInit'):
                    shutil.rmtree("DispsInit")
                f_read_obj = FcRead("fc_cart.dat")
                f_read_obj.run()
                fc_init = ForceConstant(
                    None,
                    [],
                    [],
                    0,
                    options,
                    [],
                )
                fc_init.FC =  f_read_obj.fc_mat
                os.chdir('..')
                os.chdir('..')
            else:
                # First generate displacements in internal coordinates
                eigs_init = np.eye(len(s_vec.proj.T))
                if not self.options.deriv_level:
                    indices = np.triu_indices(len(s_vec.proj.T))
                    indices = np.array(indices).T
                    print("symmetric displacements:")
                    if len(sym_sort) > 1:
                        sym_disps = []
                        for i in sym_sort:
                            print(i)
                            for j in indices:
                                if j[0] in i:
                                    if j[1] in i:
                                        sym_disps.append([j[0],j[1]])
                                        # print(j)
                        # print("sym_sort indices:")
                        # print(sym_disps)
                        # print(sym_disps[1][0])
                        # print(sym_disps[1][1])
                        # print(sym_sort)
                        # print(indices)
                        indices = sym_disps
                
                else:
                    indices = np.arange(len(eigs_init))
                if options.second_order:
                    indices = np.triu_indices(len(zmat_obj.cartesians_init.flatten()))
                    indices = np.array(indices).T
                # print("Coordinate type is:")
                # print(self.coord_type_init)
                init_disp = TransfDisp(
                    s_vec,
                    zmat_obj,
                    options.disp,
                    eigs_init,
                    True,
                    options.disp_tol,
                    TED_obj,
                    options,
                    indices,
                    deriv_level = self.options.deriv_level,
                    coord_type = self.coord_type_init
                )
                init_disp.run()
                prog_init = options.program_init
                prog_name_init = prog_init.split("@")[0]
                
                # options.calc_init = False
                if options.calc_init:
                    if os.path.exists(os.getcwd()+'/DispsInit'):
                        shutil.rmtree("DispsInit")
                    dir_obj_init = DirectoryTree(
                        prog_name_init,
                        zmat_obj,
                        init_disp,
                        options.cart_insert_init,
                        init_disp.p_disp,
                        init_disp.m_disp,
                        options,
                        indices,
                        "templateInit.dat",
                        "DispsInit",
                        deriv_level = self.options.deriv_level
                    )
                    dir_obj_init.run()
                    disp_list = []
                    for i in os.listdir(os.getcwd()):
                        disp_list.append(i)

                    if options.cluster != "sapelo":
                        v_template = VulcanTemplate(
                            options, len(disp_list), prog_name_init, prog_init
                        )
                        out = v_template.run()
                        with open("displacements.sh", "w") as file:
                            file.write(out)

                        # Submits an array, then checks if all jobs have finished every
                        # 10 seconds.
                        sub = Submit(disp_list,options)
                        sub.run()
                    else:
                        s_template = SapeloTemplate(
                            options, len(disp_list), prog_name_init, prog_init
                        )
                        out = s_template.run()
                        with open("optstep.sh", "w") as file:
                            file.write(out)
                        for z in range(0, len(disp_list)):
                            source = os.getcwd() + "/optstep.sh"
                            os.chdir("./" + str(z + 1))
                            destination = os.getcwd()
                            shutil.copy2(source, destination)
                            os.chdir("../")
                        sub = Submit(disp_list, options)
                        sub.run()

                # print("deriv_level:")
                # print(self.options.deriv_level)
                reap_obj_init = Reap(
                    options,
                    eigs_init,
                    indices,
                    options.energy_regex_init,
                    options.gradient_regex,
                    options.success_regex_init,
                    deriv_level = self.options.deriv_level
                )
                reap_obj_init.energy_regex = energy_regex
                reap_obj_init.success_regex = success_regex
                if options.calc_init:
                    reap_obj_init.run()
                    os.chdir("..")
                else:
                    os.chdir("DispsInit")
                    reap_obj_init.run()
                    os.chdir("..")

                # nate
                if not self.options.deriv_level:
                    p_array_init = reap_obj_init.p_en_array
                    m_array_init = reap_obj_init.m_en_array
                    ref_en_init = reap_obj_init.ref_en
                else:
                    cart_p_array_init = reap_obj_init.p_grad_array
                    cart_m_array_init = reap_obj_init.m_grad_array
                    p_array_init = np.zeros(np.eye(len(eigs_init)).shape)
                    m_array_init = np.zeros(np.eye(len(eigs_init)).shape)
                    ref_en_init = None
                    # Need to convert this array here from cartesians to internals using projected A-tensor
                    for i in indices:
                        grad_s_vec = SVectors(
                            zmat_obj, self.options, zmat_obj.variable_dictionary_init
                        )
                        grad_s_vec.run(init_disp.p_disp[i],False)
                        A_proj = np.dot(LA.pinv(grad_s_vec.B),TED_obj.proj)
                        p_array_init[i] = np.dot(cart_p_array_init[i].T,A_proj)
                        grad_s_vec.run(init_disp.m_disp[i],False)
                        A_proj = np.dot(LA.pinv(grad_s_vec.B),TED_obj.proj)
                        m_array_init[i] = np.dot(cart_m_array_init[i].T,A_proj)
                    

                fc_init = ForceConstant(
                    init_disp,
                    p_array_init,
                    m_array_init,
                    ref_en_init,
                    options,
                    indices,
                    deriv_level=self.options.deriv_level,
                    coord_type_init=coord_type_init
                )
                fc_init.run()
                print("Computed Force Constants:")
                print(fc_init.FC)
                if options.second_order:
                    p_array_grad = np.array([])
                    m_array_grad = np.array([])
                    for i in range(len(p_array_init)):
                        p_array_grad = np.append(p_array_grad,p_array_init[i,i])
                        m_array_grad = np.append(m_array_grad,m_array_init[i,i])
                    grad_init = ForceConstant(
                        init_disp,
                        p_array_grad,
                        m_array_grad,
                        ref_en_init,
                        options,
                        indices,
                        deriv_level=1
                    )
                    grad_init.run()
                    print("Computed Gradient:")
                    print(grad_init.FC)
                f_conv_obj = FcConv(
                    fc_init.FC,
                    s_vec,
                    zmat_obj,
                    "internal",
                    False,
                    TED_obj,
                    options.units,
                    options.second_order
                )
                f_conv_obj.N = len(fc_init.FC)
                if self.coord_type_init == "internal":
                    f_conv_obj.print_const(fc_name="fc_int_"+cma1_coord+".dat")
                    shutil.move(os.getcwd() + "/fc_int_"+cma1_coord+".dat", os.getcwd()+"/.." + self.cma1_path +"/fc_int_"+cma1_coord+".dat")
                elif self.coord_type_init == "cartesian":
                    f_conv_obj.print_const(fc_name="fc_cart.dat")
                    shutil.move(os.getcwd() + "/fc_cart.dat", os.getcwd()+"/.." + self.cma1_path +"/fc_cart.dat")
                print("Force Constants saved at:")
                print(self.cma1_path)
                if options.second_order:
                    # grad_conv_obj = FcConv(
                        # grad_init.FC,
                        # s_vec,
                        # zmat_obj,
                        # "cartesian",
                        # False,
                        # TED_obj,
                        # options.units,
                        # False
                    # )
                    fc_name = "fc_cart.grad"
                    fc_output = ""
                    g_print = grad_init.FC.copy()
                    g_print = g_print.flatten()
                    for i in range(len(g_print) // 3):
                        fc_output += "{:20.10f}".format(g_print[3 * i])
                        fc_output += "{:20.10f}".format(g_print[3 * i + 1])
                        fc_output += "{:20.10f}".format(g_print[3 * i + 2])
                        fc_output += "\n"
                    if len(g_print) % 3:
                        for i in range(len(g_print) % 3):
                            fc_output += "{:20.10f}".format(
                                g_print[3 * (len(g_print) // 3) + i]
                            )
                        fc_output += "\n"
                    with open(fc_name, "w+") as file:
                        file.write(fc_output)
                    # f_conv_obj.print_const(fc_name="fc_cart.grad")
                    shutil.move(os.getcwd() + "/fc_cart.grad", os.getcwd()+"/.." + self.cma1_path +"/fc_cart.grad")
                    print("Gradient saved at:")
                    print(self.cma1_path)
                    f_conv_obj.run(grad=grad_init.FC)
                    g_read_obj = GrRead("fc.grad")
                    g_read_obj.cart_grad = grad_init.FC

                os.chdir("..")
        
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
                options.second_order
            )
            if options.second_order:
                f_conv_obj.run(grad=g_read_obj.cart_grad)
                # f_conv_obj.F = np.dot(TED_obj.proj.T, np.dot(f_conv_obj.F, TED_obj.proj))
                # f_conv_obj.print_const(fc_name="fc_int_nat.dat")
                # raise RuntimeError
            else:
                f_conv_obj.run()
            F = f_conv_obj.F
        else:
            F = fc_init.FC
            if options.second_order:
                F = f_conv_obj.F
                F = np.dot(TED_obj.proj.T, np.dot(F, TED_obj.proj))
        
        self.options.deriv_level = 0
        
        if options.coords != "ZMAT" and not init_bool:
            F = np.dot(TED_obj.proj.T, np.dot(F, TED_obj.proj))
            if options.second_order:
                # grad_proj = np.dot(TED_obj.proj.T,grad_conv_obj.FC)
                grad_proj = np.dot(TED_obj.proj.T,f_conv_obj.v_q)
        
        if options.coords != "ZMAT":
            g_mat.G = np.dot(TED_obj.proj.T, np.dot(g_mat.G, TED_obj.proj))
        
        TED_obj.run(np.eye(TED_obj.proj.shape[1]),np.zeros(TED_obj.proj.shape[1]))
        
        # print("sym_sort:")
        # print(sym_sort)
        # if len(sym_sort) > 1:
            # Fbuff1 = np.array([])
            # Fbuff2 = {}
            # Gbuff1 = np.array([])
            # Gbuff2 = {}
            # for i in range(len(sym_sort)):
                # Fbuff1 = F.copy()
                # Fbuff1 = Fbuff1[sym_sort[i]]
                # Fbuff1 = np.array([Fbuff1[:,sym_sort[i]]])
                # Fbuff2[str(i)] = Fbuff1.copy()
                # Gbuff1 = g_mat.G.copy()
                # Gbuff1 = Gbuff1[sym_sort[i]]
                # Gbuff1 = np.array([Gbuff1[:,sym_sort[i]]])
                # Gbuff2[str(i)] = Gbuff1.copy()
            # Fbuff3 = Fbuff2[str(0)][0].copy()
            # Gbuff3 = Gbuff2[str(0)][0].copy()
            # for i in range(len(sym_sort)-1):
                # Fbuff3 = np.block([
                    # [Fbuff3,                                        np.zeros((len(Fbuff3),len(Fbuff2[str(i+1)][0])))],
                    # [np.zeros((len(Fbuff2[str(i+1)][0]),len(Fbuff3))), Fbuff2[str(i+1)][0]]
                    # ])
                # Gbuff3 = np.block([
                    # [Gbuff3,                                        np.zeros((len(Gbuff3),len(Gbuff2[str(i+1)][0])))],
                    # [np.zeros((len(Gbuff2[str(i+1)][0]),len(Gbuff3))), Gbuff2[str(i+1)][0]]
                    # ])
            # F = Fbuff3
            # g_mat.G = Gbuff3
        if len(sym_sort) > 1:
            print("Initial Force Constants:")
            print(flat_sym_sort)
            F_sym = F[flat_sym_sort].copy()
            F_sym = F_sym[:,flat_sym_sort]
            print(F_sym)
            
            print("Initial G-Matrix:")
            g_sym = g_mat.G[flat_sym_sort].copy()
            g_sym = g_sym[:,flat_sym_sort]
            g_sym[np.abs(g_sym) < 1e-9] = 0
            print(sym_sort)
            print(g_sym)

        
        print("Initial Force Constants:")
        print(F.shape)
        print(F)
        
        print("Initial G-Matrix:")
        g_mat.G[np.abs(g_mat.G) < 1e-9] = 0
        print(g_mat.G)

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
        
        print("TED for sym purposes: ")
        print(init_GF.ted.TED)
        
        self.ref_init = init_GF.freq
        if len(sym_sort):
            ref_TED_init = init_GF.ted.TED
            self.irreps_init = []
            # sym_freqs_ref_init = []
            for irrep in sym_sort:
                print("Sym natty indices:")
                print(irrep)
                irrep_modes = []
                for i in range(len(ref_TED_init.T)):
                    Sum = 0
                    for j in irrep:
                        # print(j)
                        Sum += ref_TED_init.T[i,j]
                    print(Sum)
                    if Sum > 98.:
                        irrep_modes.append(i)
                print(irrep_modes)
                if len(irrep_modes) != len(irrep):
                    print("Something's wrong with the irrep symmetry sorter:")
                    raise RuntimeError
                self.irreps_init.append(irrep_modes)
                # sym_freqs_ref_init.append(irrep_modes)
            print("Modes sorted by sym:")
            print(self.irreps_init)
            print("Sym_sort for reference:")
            print(sym_sort)

            sym_freqs_ref_init = copy.deepcopy(self.irreps_init)
            del_list = []
            for i in range(len(self.irreps_init)):
                if len(self.irreps_init[i]) > 1:
                    for j in range(len(self.irreps_init[i])):
                        index = self.irreps_init[i][j]
                        sym_freqs_ref_init[i][j] = self.ref_init[index].copy()
                    sym_freqs_ref_init[i].reverse()
                else:
                    del_list.append(i)
            print(sym_freqs_ref_init)
            print(self.irreps_init)
            del_list.reverse()
            for i in del_list:
                del sym_freqs_ref_init[i]
            print("Initial Freqs:")
            print(init_GF.freq)
            print("Sorted and trimmed freqs:")
            print(sym_freqs_ref_init)
            flat_sym_freq_init = [
                x
                for xs in sym_freqs_ref_init
                for x in xs
            ]
            self.ref_init = np.array(flat_sym_freq_init)
            print(self.ref_init)
            # raise RuntimeError


        # Now for the TED check.
        G = np.dot(np.dot(LA.inv(init_GF.L), g_mat.G), LA.inv(init_GF.L).T)
        G[np.abs(G) < options.tol] = 0
        F = np.dot(np.dot(init_GF.L.T, F), init_GF.L)
        F[np.abs(F) < options.tol] = 0
        # if self.options.second_order:
            # grad_n = np.dot(init_GF.L.T, grad_proj)
            # print("Normal Mode Gradients:")
            # for i in range(len(grad_n)):
                # print(str(i+1) + ": " + str(grad_n[i]))
        
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
       
        # eigs = len(TED_GF.S)
        # print('eigs')
        # print(eigs)
        # self.eigs = eigs


 
        proj_tol = 1.0e-3
        eig_inv = inv(init_GF.L)  # (Normal modes (Q) x Sym internals (S) )
        for i in range(len(eig_inv)):
            eig_inv[i] = eig_inv[i] / LA.norm(eig_inv[i])
            eig_inv[i][
                np.abs(eig_inv[i]) < np.max(np.abs(eig_inv[i])) * proj_tol
            ] = 0
        
        # Now run the TZ force constant transformation
        print(os.getcwd())
        zmat_obj2 = Zmat(options)
        zmat_obj2.run(zmat_name="zmat2")
        
        print(options.man_proj)
        options.man_proj = True
        
        
        s_vec = SVectors(
            zmat_obj2, options, zmat_obj2.variable_dictionary_init
        )
        s_vec.run(zmat_obj2.cartesians_init, True, proj=TED_obj.proj)
                
        TED_obj = TED(s_vec.proj, zmat_obj2)
                
        g_mat = GMatrix(zmat_obj2, s_vec, options)
        g_mat.run()
        
        # if len(sym_sort) > 1:
            # flat_sym_sort = np.array([])
            # for i in range(len(sym_sort)):
                # flat_sym_sort = np.append(flat_sym_sort,sym_sort[i])
            # flat_sym_sort = flat_sym_sort.astype(int)
        
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
                False,
            )
            f_conv_obj.run()
            F = f_conv_obj.F
        else:
            F = fc_init.FC
        
        
        
        # redundant basis 
        G = np.dot(np.dot(TED_obj.proj.T,G),TED_obj.proj)
        # if len(sym_sort) > 1:
            # G = G[flat_sym_sort]
            # G = G[:,flat_sym_sort]
        print("Giraffe G")
        G[np.abs(G) < 1.0e-9] = 0 
        print(G)
        if len(sym_sort) > 1:
            print(g_sym)
        G = np.dot(np.dot(eig_inv, G), eig_inv.T)
        # G[np.abs(G) < options.tol] = 0
        # print("Nat F:")
        # F[np.abs(F) < 1.0e-5] = 0 
        # print(F)
        # Conversion to aJ/Ang
        F_aJ = F.copy()
        F_aJ *= 4.3597447222071
        F_aJ /= 0.529177210903
        F = np.dot(np.dot(TED_obj.proj.T,F),TED_obj.proj)
        F_aJ = np.dot(np.dot(TED_obj.proj.T,F_aJ),TED_obj.proj)
        # if len(sym_sort) > 1:
            # F = F[flat_sym_sort]
            # F = F[:,flat_sym_sort]
            # F_aJ = F_aJ[flat_sym_sort]
            # F_aJ = F_aJ[:,flat_sym_sort]
        print("Giraffe F")
        F[np.abs(F) < 1.0e-5] = 0 
        print(F)
        if len(sym_sort) > 1:
            print(sym_sort)
            print(F_sym)
            F_sym = F[flat_sym_sort]
            F_sym = F_sym[:,flat_sym_sort]
        # print("aJ F")
        # F_aJ[np.abs(F) < 1.0e-5] = 0 
        # print(F_aJ)
        # if len(sym_sort) > 1:
            # print("aJ/A Sym Force Constants:")
            # print(flat_sym_sort+1)
            # F_aJ = F_aJ[flat_sym_sort]
            # F_aJ = F_aJ[:,flat_sym_sort]
            # print(F_aJ)
        F = np.dot(np.dot(inv(eig_inv).T, F), inv(eig_inv))
        # F[np.abs(F) < options.tol] = 0
         
        full_GF = GFMethod(
            G,
            F,
            options.tol,
            options.proj_tol,
            zmat_obj2,
            TED_obj,
            False
        )
        full_GF.run()
        self.ted = full_GF.ted.TED      # TED matrix
        
        # Print Full TED here in projected basis

        print("////////////////////////////////////////////")
        print("//{:^40s}//".format(" Full Hessian TED"))
        print("////////////////////////////////////////////")
        TED_obj.run(np.dot(init_GF.L, full_GF.L), full_GF.freq, rect_print=False)
        
        self.reference_freq = full_GF.freq 
        if len(sym_sort):
            ref_TED = TED_obj.TED
            self.irreps_ref = []
            # sym_freqs_ref_init = []
            for irrep in sym_sort:
                print("Sym natty indices:")
                print(irrep)
                irrep_modes = []
                for i in range(len(ref_TED.T)):
                    print("Mode index:")
                    print(i)
                    Sum = 0
                    for j in irrep:
                        # print(j)
                        Sum += ref_TED.T[i,j]
                    print(Sum)
                    if Sum > 98.:
                        irrep_modes.append(i)
                if len(irrep_modes) != len(irrep):
                    print("Something's wrong with the irrep symmetry sorter:")
                    raise RuntimeError
                self.irreps_ref.append(irrep_modes)
                # sym_freqs_ref_init.append(irrep_modes)
            print("Ref modes sorted by sym:")
            print(self.irreps_ref)
            print("Sym_sort for reference:")
            print(sym_sort)

            sym_freqs_ref = copy.deepcopy(self.irreps_ref)
            del_list = []
            for i in range(len(self.irreps_ref)):
                if len(self.irreps_ref[i]) > 1:
                    for j in range(len(self.irreps_ref[i])):
                        index = self.irreps_ref[i][j]
                        sym_freqs_ref[i][j] = self.reference_freq[index].copy()
                    sym_freqs_ref[i].reverse()
                else:
                    del_list.append(i)
            print(sym_freqs_ref)
            print(self.irreps_ref)
            del_list.reverse()
            for i in del_list:
                del sym_freqs_ref[i]
            print("Initial Freqs:")
            print(init_GF.freq)
            print("Sorted and trimmed freqs:")
            print(sym_freqs_ref)
            flat_sym_freq = [
                x
                for xs in sym_freqs_ref
                for x in xs
            ]
            self.reference_freq = np.array(flat_sym_freq)
            print(self.reference_freq)
            # raise RuntimeError
        
        m = 2 
        var = 0.95 
        
        # def checkted(ted):
            # temps = []
            # for i in range(0,np.shape(ted)[0]):
                # ted_slice = ted[:,i] 
                # temp = copy.copy(ted_slice)
                # for j in range(0,m):
                    # largest = np.argmax(temp)
                    # if temp[largest] < 0.9:
                        # print('not big enough')
                    # print('largest')
                    # print(largest,temp[largest])
                    # temps.append([i,largest])
                    # temp[largest] = 0
                    # print('another')
                    # print(temp)
            # return temps
        # print('is this the ted im looking for?')
        # ted_breakdown = init_GF.ted_breakdown
        # print(ted_breakdown)  
        
        # temps = checkted(ted_breakdown) 
        
        # print('temps')
        # print(temps)
        
        self.reference_TED = TED_obj.TED
        self.exclude = np.array([])
        # print(sym_sort)
        # print(self.lone_wolves)
        # print(self.reference_TED)
        ref_TED = self.reference_TED
        for wolf in self.lone_wolves:
            for i in range(len(ref_TED.T)):
                for j in range(len(ref_TED.T[i])):
                    if j == wolf and ref_TED.T[i,j] > 99.0:
                        # print(i,j)
                        # print(ref_TED.T[i,j])
                        self.exclude = np.append(self.exclude,i)
        self.exclude = self.exclude.astype(int)

        # irreps = []
        # for irrep in sym_sort:
            # print("Irrep indices:")
            # print(irrep)
            # irrep_modes = []
            # for i in range(len(ref_TED.T)):
                # # print("Mode index:")
                # # print(i)
                # Sum = 0
                # for j in irrep:
                    # # print(j)
                    # Sum += ref_TED.T[i,j]
                # # print(Sum)
                # if Sum > 99.:
                    # irrep_modes.append(i)
            # print(irrep_modes)
            # if len(irrep_modes) != len(irrep):
                # print("Something's wrong with the irrep symmetry sorter:")
                # raise RuntimeError
            # irreps.append(irrep_modes)
        # print("Modes sorted by sym:")
        # print(irreps)
        # print("Sym_sort for reference:")
        # print(sym_sort)


        # print(self.exclude)
        # raise RuntimeError
        if options.coords == 'Redundants':
            L_B = full_GF.L
        elif options.coords == 'Custom':
            L_A = full_GF.L
         
        
        def n_largest(n, FC):
            indexes = []
            upper_triang = abs(np.triu(FC,n))
            for i in range(0,n):
                fc_cma2 = np.where(upper_triang == upper_triang.max())
                index = [fc_cma2[0][0], fc_cma2[1][0]]
                indexes.append(index)
                upper_triang[index[0],index[1]] = 0
            print(indexes)
            return indexes
        
        np.set_printoptions(edgeitems=60,linewidth=10000)
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
        Fdiag = copy.copy(np.diag(np.diag(F)))
        
        print("Diagonal Force constant matrix in lower level normal mode basis:")
        print(Fdiag)
        diag_GF = GFMethod(
            G,
            Fdiag,
            options.tol,
            options.proj_tol,
            zmat_obj2,
            TED_obj,
            False
        )
        
        diag_GF.run()

        self.Freq_CMA0 = diag_GF.freq

        # if options.coords == 'Redundant':
            # self.Freq_redundant = diag_GF.freq
        # elif options.coords == 'Custom':
            # self.Freq_custom = diag_GF.freq
        # elif options.coords == 'ZMAT' :
            # self.Freq_zmat = diag_GF.freq
        # else:
            # pass
        
        
        # Print Diagonal TED here in projected basis

        print("////////////////////////////////////////////")
        print("//{:^40s}//".format(" CMA-0 TED"))
        print("////////////////////////////////////////////")
        TED_obj.run(np.dot(init_GF.L, diag_GF.L), diag_GF.freq, rect_print=False)
        
        if len(sym_sort):
            ref_TED = TED_obj.TED
            self.irreps_CMA0 = []
            for irrep in sym_sort:
                print("Sym natty indices:")
                print(irrep)
                irrep_modes = []
                for i in range(len(ref_TED.T)):
                    print("Mode index:")
                    print(i)
                    Sum = 0
                    for j in irrep:
                        # print(j)
                        Sum += ref_TED.T[i,j]
                    print(Sum)
                    if Sum > 98.:
                        irrep_modes.append(i)
                if len(irrep_modes) != len(irrep):
                    print("Something's wrong with the irrep symmetry sorter:")
                    raise RuntimeError
                self.irreps_CMA0.append(irrep_modes)
            print("Ref modes sorted by sym:")
            print(self.irreps_CMA0)
            print("Sym_sort for reference:")
            print(sym_sort)

            sym_freqs_CMA0 = copy.deepcopy(self.irreps_CMA0)
            del_list = []
            for i in range(len(self.irreps_CMA0)):
                if len(self.irreps_CMA0[i]) > 1:
                    for j in range(len(self.irreps_CMA0[i])):
                        index = self.irreps_CMA0[i][j]
                        sym_freqs_CMA0[i][j] = self.Freq_CMA0[index].copy()
                    sym_freqs_CMA0[i].reverse()
                else:
                    del_list.append(i)
            print(sym_freqs_CMA0)
            print(self.irreps_CMA0)
            del_list.reverse()
            for i in del_list:
                del sym_freqs_CMA0[i]
            print("Initial Freqs:")
            print(init_GF.freq)
            print("Sorted and trimmed freqs:")
            print(sym_freqs_CMA0)
            flat_sym_freq = [
                x
                for xs in sym_freqs_CMA0
                for x in xs
            ]
            self.Freq_CMA0 = np.array(flat_sym_freq)
            print(self.Freq_CMA0)
            # raise RuntimeError


        if self.options.n_cma2 > 0:
            self.Freq_cma2 = np.array([])
            self.eta_num = np.array([])
            self.eta_denom = np.array([])
            self.total_off_diags = np.array([])
            for xi_tol_i in xi_tol:
                if self.options.off_diag:
                    print('this is the option')
                    print(self.options.off_diag)
                    algo = Algorithm(eigs, None, options)
                    algo.run()
                    print('algo indices')
                    print(algo.indices)
                    temp = np.zeros((eigs,eigs))  
                    print('temp')
                    print(temp) 
                    for z, extra in enumerate(algo.indices):
                        element = F[extra[0], extra[1]] 
                        temp[extra[0], extra[1]] = element
                        temp[extra[1], extra[0]] = element
                    print('temp')
                    print(temp)
                else:
                    # extras = n_largest(self.options.n_cma2, np.abs(copy.copy(F)))
                    #91_ccsd extras = [[4,5],[10,11],[14,16]] 
                    #72_ccsd extras = [[4,5]]
                    #90_ccsd extras = [[15,17]]
                    #2.16_ccsd extras = [[7,9],[8,9]]
                    #10_ccsd extras = [[12,15]]
                    #82_ccsd extras = [[14,15]]
                    #70_ccsd extras = [[10,11],[0,1]]
                    #85_ccsd extras = [[1,4],[3,5]]
                    extras = []
                    # extras = [[16,17],[17,18],[14,17],[16,18],[14,16]]
                    # extras = [[4,5],[10,11],[14,16],[16,17],[16,18]]
                    extras = [[7,8]]
                    # print("Inter off diag indices:")
                    # for i in range(len(Fdiag)):
                        # for j in range(len(Fdiag)-i-1):
                            # # print([i,j+i+1])
                            # extras.append([i,j+i+1])

                    # extras = [[0,1],[0,2],[0,3],[0,4],[0,5],[0,6],[0,7],[0,8],[1,2],[1,3],[1,4],[1,5],[1,6],[1,7],[1,8],[2,3],[2,4],[2,5],[2,6],[2,7],[2,8],[3,4],[3,5],[3,6],[3,7],[3,8],[4,5],[4,6],[4,7],[4,8],[5,6],[5,7],[5,8],[6,7],[6,8],[7,8]]
                    # extras = [[0,1],[0,2],[0,3],[0,4],[0,5],[1,2],[1,3],[1,4],[1,5],[2,3],[2,4],[2,5],[3,4],[3,5],[4,5]]
                    # print(np.array(extras).shape)
                    # extras = [[0,1],[0,2],[0,3],[0,4],[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]]
                    #extras = [[17,19]]
                    #extras = [[0,2]]
                    print('extras')
                    print(extras)
                    # raise RuntimeError
                    temp = copy.copy(Fdiag)
                    if len(self.options.other_F_matrix) and os.path.exists(os.getcwd() + "/inter_fc.dat"):
                        f_read_obj_inter = FcRead("inter_fc.dat")
                        f_read_obj_inter.run()
                        F_inter = f_read_obj_inter.fc_mat
                        F_inter = np.dot(np.dot(inv(eig_inv).T, F_inter), inv(eig_inv))
                        print("F_inter:")
                        print(F_inter)
                        print("F_A:")
                        print(F)
                        xi = F_inter * 0
                        extras = []
                        if len(sym_sort) > 1:
                            self.total_off_diags_buff = 0
                            print(self.irreps_init)
                            for irrep in self.irreps_init:
                                print(irrep)
                                if len(irrep) > 1:
                                    for i in range(len(irrep)):
                                        for j in range(i):
                                            if i != j:
                                                print(i,j)
                                                a = irrep[i]
                                                b = irrep[j]
                                                buff = np.abs(F_inter[a,b])
                                                xi[a,b] = buff / np.sqrt(np.abs(F_inter[a,a])*np.abs(F_inter[b,b]))
                                                if xi[a,b] > xi_tol_i:
                                                    print(a,b)
                                                    print(xi[a,b])
                                                    extras.append([a,b])
                                    self.total_off_diags_buff += (len(irrep)**2 - len(irrep))/2
                        else:
                            self.total_off_diags_buff = (len(xi)**2 - len(xi))/2
                            for i in range(len(xi)):
                                for j in range(i+1):
                                    if i != j:
                                        # print(i,j)
                                        buff = np.abs(F_inter[i,j])
                                        xi[i,j] = buff / np.sqrt(np.abs(F_inter[i,i])*np.abs(F_inter[j,j]))
                                        if xi[i,j] > xi_tol_i:
                                            print(i,j)
                                            print(xi[i,j])
                                            extras.append([i,j])
                        
                        self.cma_off_diags = len(extras)
                        # self.off_diags = np.append(self.off_diags,len(extras))
                        self.total_off_diags = np.append(self.total_off_diags,self.total_off_diags_buff)
                        # self.perc_off_diags = (self.cma_off_diags / self.total_off_diags) * 100
                        self.eta_num = np.append(self.eta_num,self.cma_off_diags*1.0)
                        self.eta_denom = np.append(self.eta_denom,len(xi)*1.0)
                        # raise RuntimeError
                        for z, extra in enumerate(extras):
                            element = F[extra[0], extra[1]] 
                            temp[extra[0], extra[1]] = element
                            temp[extra[1], extra[0]] = element
                        
                        
                    if len(self.options.other_F_matrix) and os.path.exists(os.getcwd() + "/inter_fc.dat"):
                        pass
                        # print(os.getcwd())
                        # f_read_obj_inter = FcRead("inter_fc.dat")
                        # f_read_obj_inter.run()
                        # F_inter = f_read_obj_inter.fc_mat
                        # F_inter = np.dot(np.dot(inv(eig_inv).T, F_inter), inv(eig_inv))
                        # print("F_inter:")
                        # print(F_inter)
                        # print("F_A:")
                        # print(F)
                        # # raise RuntimeError
                        # for z, extra in enumerate(extras):
                            # element = F_inter[extra[0], extra[1]] 
                            # temp[extra[0], extra[1]] = element
                            # temp[extra[1], extra[0]] = element
                    else:
                        for z, extra in enumerate(extras):
                            element = F[extra[0], extra[1]] 
                            temp[extra[0], extra[1]] = element
                            temp[extra[1], extra[0]] = element
                    print('CMA2 FC matrix')
                    print(temp) 
                #if options.coords == 'Redundant':
                #    #F[index] = self.F_redundant[index]     
                #if options.coords == 'Custom':
                #    extras = n_largest(2, np.abs(copy.copy(self.F_custom)))
                #    #F[index] = self.F_custom[index]     
                #elif options.coords == 'ZMAT' :
                #    extras = n_largest(2, np.abs(copy.copy(self.F_zmat)))
                #    #F[index] = self.F_zmat[index]     
                #else:
                #    pass
                print('Time for some off-diags')
                cma2_GF = GFMethod(
                    G,
                    temp,
                    options.tol,
                    options.proj_tol,
                    zmat_obj2,
                    TED_obj,
                    False
                )
                cma2_GF.run()
                # print('CMA2 including ' + str(z + 1) + ' off-diagonal bands/elements for ' + str(options.coords) + ' coordinates')
                cma2_Freq = cma2_GF.freq.copy()
                # cma2_Freq = np.delete(cma2_Freq,self.exclude)

                
                print("////////////////////////////////////////////")
                print("//{:^40s}//".format(" CMA-2 TED"))
                print("////////////////////////////////////////////")
                TED_obj.run(np.dot(init_GF.L, cma2_GF.L), cma2_GF.freq, rect_print=False)
                
                if len(sym_sort):
                    ref_TED = TED_obj.TED
                    self.irreps_CMA2 = []
                    for irrep in sym_sort:
                        print("Sym natty indices:")
                        print(irrep)
                        irrep_modes = []
                        for i in range(len(ref_TED.T)):
                            print("Mode index:")
                            print(i)
                            Sum = 0
                            for j in irrep:
                                # print(j)
                                Sum += ref_TED.T[i,j]
                            print(Sum)
                            if Sum > 98.:
                                irrep_modes.append(i)
                        if len(irrep_modes) != len(irrep):
                            print("Something's wrong with the irrep symmetry sorter:")
                            raise RuntimeError
                        self.irreps_CMA2.append(irrep_modes)
                    print("Ref modes sorted by sym:")
                    print(self.irreps_CMA2)
                    print("Sym_sort for reference:")
                    print(sym_sort)

                    sym_freqs_CMA2 = copy.deepcopy(self.irreps_CMA2)
                    del_list = []
                    for i in range(len(self.irreps_CMA2)):
                        if len(self.irreps_CMA2[i]) > 1:
                            for j in range(len(self.irreps_CMA2[i])):
                                index = self.irreps_CMA2[i][j]
                                sym_freqs_CMA2[i][j] = cma2_Freq[index].copy()
                            sym_freqs_CMA2[i].reverse()
                        else:
                            del_list.append(i)
                    print(sym_freqs_CMA2)
                    print(self.irreps_CMA2)
                    del_list.reverse()
                    for i in del_list:
                        del sym_freqs_CMA2[i]
                    print("Sorted and trimmed freqs:")
                    print(sym_freqs_CMA2)
                    flat_sym_freq = [
                        x
                        for xs in sym_freqs_CMA2
                        for x in xs
                    ]
                    cma2_Freq = np.array(flat_sym_freq)
                    print(cma2_Freq)

            self.Freq_cma2 = np.append(self.Freq_cma2,cma2_Freq,axis=0)
            self.Freq_cma2 = np.reshape(self.Freq_cma2,(len(xi_tol),-1))
            # print("CMA2 Freqs:")
            # print(cma2_GF.freq) 
        
        def n_largest(n, FC):
            indexes = []
            upper_triang = abs(np.triu(FC,1))
            #print('this is the upper triang')
            #print(upper_triang)
            length = len(upper_triang)
            for i in range(0,n):
                index = np.argmax(upper_triang)
                if index > length:
                    two_d = [index // length, index % length]
                else:
                    two_d = [0,index]
                indexes.append(two_d)
                
                upper_triang[two_d[0],two_d[1]] = 0
            return indexes

