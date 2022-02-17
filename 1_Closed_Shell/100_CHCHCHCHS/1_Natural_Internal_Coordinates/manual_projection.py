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

        a, b = np.cos(144*np.pi/180), np.cos(72*np.pi/180)

        HA_str1 = np.eye(1)
       
        HA_str2 = normalize(np.array([
        [1, 1],
        [1,-1],
        ]).T)
       
        HA_str3 = normalize(np.array([
        [1, 1],
        [1,-1],
        ]).T)

        CH_str1 = normalize(np.array([
        [1, 1],
        [1,-1],
        ]).T)
       
        CH_str2 = normalize(np.array([
        [1, 1],
        [1,-1],
        ]).T)

        HA_ang = normalize(np.array([
        [1,   a,   a,   b,   b],
        [0, a-b, b-a, 1-a, a-1],
        ]).T)

        CH_ang1 = normalize(np.array([
        [1,-1, 1,-1],
        [1,-1,-1, 1],
        ]).T)

        CH_ang2 = normalize(np.array([
        [1,-1, 1,-1],
        [1,-1,-1, 1],
        ]).T)

        tor = normalize(np.array([
        [1,   b,   b,   a,   a],
        [0, 1-a, a-1, a-b, b-a],
        ]).T)

        oop2 = normalize(np.array([
        [1, 1],
        [1,-1],
        ]).T)

        oop3 = normalize(np.array([
        [1, 1],
        [1,-1],
        ]).T)

        Proj = block_diag(HA_str1,HA_str2,HA_str3,CH_str1,CH_str2,HA_ang,CH_ang1,CH_ang2,tor,oop2,oop3)

        self.Proj = Proj

def normalize(mat):
    return 1/norm(mat,axis=0)*mat

if __name__=="__main__":
    np.set_printoptions(linewidth=400, precision=2,threshold=100000)
    p = Projection([])
    p.run()
    print(p.Proj)

