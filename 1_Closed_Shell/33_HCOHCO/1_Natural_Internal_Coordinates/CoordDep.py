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

options_kwargs = {
    "queue": "gen4.q,gen6.q,gen5.q",
    "program": "molpro@2010.1.67+mpi",
    "energy_regex": r"\(T\) energy\s+(\-\d+\.\d+)",
    'energy_regex' : r"\s*\!CCSD\(T\) total energy\s+(-\d+\.\d+)",
    "cart_insert": 9,
    "man_proj" : True,
    "calc" : False,
    "coords": "Custom",
    "success_regex": r"Variable memory released",
}
options_obj = Options(**options_kwargs)



options = options_obj

rootdir = os.getcwd()
zmat_obj = Zmat(options)
zmat_obj.run()
        
# Build manual projection matrix here.
np.set_printoptions(edgeitems=60,linewidth=1000,precision=4)
""" 
    Thus begins the code to modify
"""

CC_mat = np.eye(1)

CO_mat = (1/np.sqrt(2)) * np.array([
    [1,1],
    [1,-1]
])

angles_mat = np.array([
    [0.5,0.5,0.5,0.5],
    [0.5,-0.5,-0.5,0.5],
    [0.5,0.5,-0.5,-0.5],
    [0.5,-0.5,0.5,-0.5]
])

oop_mat = np.array([
    [1/np.sqrt(2),1/np.sqrt(2)],
    [1/np.sqrt(2),-1/np.sqrt(2)]
])

tor_mat = np.array([
    [1/np.sqrt(2)],
    [1/np.sqrt(2)]
])

Proj = block_diag(CC_mat,CO_mat,CO_mat,angles_mat,tor_mat,oop_mat)
""" 
    Thus ends the code to modify
"""
print(Proj)
# print(Proj.shape)


# Compute the initial s-vectors
s_vec = SVectors(
    zmat_obj, options, zmat_obj.variable_dictionary_init
)
s_vec.run(zmat_obj.cartesians_init, True, proj=Proj)
        
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


# Now run the TZ force constant transformation

zmat_obj2 = Zmat(options)
zmat_obj2.run(zmat_name="zmat2")


s_vec = SVectors(
    zmat_obj2, options, zmat_obj2.variable_dictionary_init
)
s_vec.run(zmat_obj2.cartesians_init, True, proj=Proj)
        
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

G = np.dot(np.dot(TED_obj.proj.T,G),TED_obj.proj)
G = np.dot(np.dot(eig_inv, G), eig_inv.T)
# G[np.abs(G) < options.tol] = 0
F = np.dot(np.dot(TED_obj.proj.T,F),TED_obj.proj)
F = np.dot(np.dot(inv(eig_inv).T, F), inv(eig_inv))
F[np.abs(F) < options.tol] = 0


print("Full Force constant matrix in lower level normal mode basis:")
print(F)
F = np.diag(np.diag(F))

print("Diagonal Force constant matrix in lower level normal mode basis:")
print(F)

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
