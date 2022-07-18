import numpy as np
from numpy.linalg import norm
from scipy.linalg import block_diag

class Projection(object):
    """
    This class is used to specify the manual projection matrix
    for CMA. It is stored as an object and is only needed when
    self.options.man_proj = True.
    """

    def __init__(self,  options):

        self.options = options

    def run(self):

        HA_str = normalize(np.array([
        [ 1, 1, 1],
        [ 2,-1,-1],
        [ 0, 1,-1],
        ]).T)
       
        NH_str = normalize(np.array([
        [1, 1, 1, 1, 1, 1, 1, 1, 1],
        [2, 2, 2,-1,-1,-1,-1,-1,-1],
        [0, 0, 0, 1, 1, 1,-1,-1,-1],
        [2,-1,-1, 2,-1,-1, 2,-1,-1],
        [4,-2,-2,-2, 1, 1,-2, 1, 1],
        [0, 0, 0, 2,-1,-1,-2, 1, 1],
        [0, 1,-1, 0, 1,-1, 0, 1,-1],
        [0, 2,-2, 0,-1, 1, 0,-1, 1],
        [0, 0, 0, 0, 1,-1, 0,-1, 1],
        ]).T)

        HA_ang = normalize(np.array([
        [2,-1,-1],
        [0, 1,-1],
        ]).T)

        NH_ang = normalize(np.array([
        [1, 1, 1,-1,-1,-1, 1, 1, 1,-1,-1,-1, 1, 1, 1,-1,-1,-1],
        [2, 2, 2,-2,-2,-2,-1,-1,-1, 1, 1, 1,-1,-1,-1, 1, 1, 1],
        [0, 0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1,-1,-1,-1, 1, 1, 1],
        [2,-1,-1, 0, 0, 0, 2,-1,-1, 0, 0, 0, 2,-1,-1, 0, 0, 0],
        [4,-2,-2, 0, 0, 0,-2, 1, 1, 0, 0, 0,-2, 1, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 2,-1,-1, 0, 0, 0,-2, 1, 1, 0, 0, 0],
        [0, 1,-1, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1, 0, 0, 0],
        [0, 2,-2, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0,-1, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0,-1, 1, 0, 0, 0],
        [0, 0, 0, 2,-1,-1, 0, 0, 0, 2,-1,-1, 0, 0, 0, 2,-1,-1],
        [0, 0, 0, 4,-2,-2, 0, 0, 0,-2, 1, 1, 0, 0, 0,-2, 1, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-1,-1, 0, 0, 0,-2, 1, 1],
        [0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1],
        [0, 0, 0, 0, 2,-2, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0,-1, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0,-1, 1],
        ]).T)

        tor = normalize(np.array([
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        [2, 2, 2, 2, 2, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
        [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1],
        ]).T) 

        buff = np.array([
        [0, 0, 0],
        ]).T 
        oop = normalize(np.array([
        [ 0.,
          0.,
          0.,
          0.,
          0.,
          0.,
          0.,
          0.,
          0.,
          0.,
          0.,
          0.,
          1.,
          1.,
          1.,
          0.,
          0.,
          0.,
          0.,
          0.,
          0.,
          0.,
          0.,
          0.,
          0.,
          0.,
          0.,
          0.,
          0.,
          0.,
          0.,
          0.,
          0.,
          0.,
          0.,
          0.,
          0.,
          0.,
          0.,
          0.,
          0.,
          0.,
          0.,
          0.,
          0.,
          0.,
          0.,
          0.,
          0.,
          0.,
          0.,
         -1.,
         -1.,
         -1.]
        ]).T) 

        Proj = block_diag(HA_str,NH_str,HA_ang,NH_ang,tor,buff)
        # Proj = block_diag(HA_str,NH_str,HA_ang,NH_ang,tor)
        projBuff = Proj.copy()
        projBuff = np.append(projBuff[:,:14],oop,axis=1)
        Proj = np.append(projBuff,Proj[:,14:-1],axis=1)
        self.sym_sort = np.array([
            [0,
            3,
            6,
            14,
            15,
            18,
            24],
            [9,
            21,
            27,
            30],
            [1,
            2,
            4,
            5,
            7,
            8,
            10,
            11,
            12,
            13,
            16,
            17,
            19,
            20,
            22,
            23,
            25,
            26,
            28,
            29,
            31,
            32]
            ],dtype=object)

        self.Proj = Proj
        np.set_printoptions(edgeitems=60,linewidth=1000)
        # print(Proj.T)
        # raise RuntimeError

def normalize(mat):
    return 1/norm(mat,axis=0)*mat

if __name__=="__main__":
    np.set_printoptions(linewidth=400, precision=2,threshold=100000)
    p = Projection([])
    p.run()
    print(p.Proj)

