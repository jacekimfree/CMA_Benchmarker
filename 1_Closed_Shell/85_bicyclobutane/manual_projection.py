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

        # HA_str = normalize(np.array([
        # [2, 1, 1, 1, 1],
        # [0, 1, 1, -1,-1],
        # [4,-1,-1, -1, -1],
        # [0, 1,-1, 1,-1],
        # [0, 1,-1,-1, 1],
        # ]).T)
        HA_str = normalize(np.array([
        [ 1, 1, 1, 1],
        [ 1, 1, -1, -1],
        [ 1,-1, 1,-1],
        [ 1,-1,-1, 1],
        ]).T)
       
        CH_str1 = normalize(np.array([
        [1, 1],
        [1,-1],
        ]).T)
       
        CH_str2 = normalize(np.array([
        [1, 1, 1, 1],
        [1, 1,-1,-1],
        [1,-1, 1,-1],
        [1,-1,-1, 1],
        ]).T)
       
        CH_ang1 = normalize(np.array([
        [1,-1, 1,-1],
        [1,-1,-1, 1],
        ]).T)

        CH_ang2 = normalize(np.array([
        [4,-1,-1,-1,-1, 4,-1,-1,-1,-1],
        [4,-1,-1,-1,-1,-4, 1, 1, 1, 1],
        [0, 1, 1,-1,-1, 0, 1, 1,-1,-1],
        [0, 1, 1,-1,-1, 0,-1,-1, 1, 1],
        [0, 1,-1, 1,-1, 0, 1,-1, 1,-1],
        [0, 1,-1, 1,-1, 0,-1, 1,-1, 1],
        [0, 1,-1,-1, 1, 0, 1,-1,-1, 1],
        [0, 1,-1,-1, 1, 0,-1, 1, 1,-1],
        ]).T)

        tor = np.eye(1)

        oop = normalize(np.array([
        [1, 1],
        [1,-1],
        ]).T)

        # Proj = block_diag(HA_str,CH_str1,CH_str2,CH_ang1,CH_ang2,tor,oop)
        Proj = block_diag(np.eye(1),HA_str,CH_str1,CH_str2,CH_ang1,CH_ang2,tor,oop)

        self.Proj = Proj

def normalize(mat):
    return 1/norm(mat,axis=0)*mat

if __name__=="__main__":
    np.set_printoptions(linewidth=400, precision=2,threshold=100000)
    p = Projection([])
    p.run()
    print(p.Proj)

