import numpy as np
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

        str1 = np.array([
            [0.5, 0.5, 0.5, 0.5],
            [0.5, 0.5,-0.5,-0.5],
            [1/np.sqrt(2),-1/np.sqrt(2), 0, 0],
            [0, 0, 1/np.sqrt(2),-1/np.sqrt(2)]
        ]).T
        
        ang1 = 1/np.sqrt(2) * np.array([
            [1, 1],
            [1,-1]
        ])
        
        ang2 = 0.5 * np.array([
            [ 1, 1, 1],
            [-1, 1,-1],
            [ 1,-1,-1],
            [-1,-1, 1],
        ])
        
        Proj = block_diag(str1, ang1, ang2)
        
        
        self.Proj = Proj                     
