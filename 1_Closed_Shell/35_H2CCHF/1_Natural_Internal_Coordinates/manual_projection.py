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

        CC_mat = np.eye(1)
        
        CX_mat = np.eye(2)
        
        CH_mat = np.array([
             [1/np.sqrt(2), 1/np.sqrt(2)],
             [1/np.sqrt(2),-1/np.sqrt(2)]
        ])
        
        CX_ang_mat = np.array([
            [-1/np.sqrt(6), 1/np.sqrt(2)],
            [ 2/np.sqrt(6), 0           ],
            [-1/np.sqrt(6),-1/np.sqrt(2)]
        ])
        
        CH_ang_mat = np.array([
            [ 2/np.sqrt(6), 0           ],
            [-1/np.sqrt(6), 1/np.sqrt(2)],
            [-1/np.sqrt(6),-1/np.sqrt(2)]
        ])
        
        tor_mat = 1/np.sqrt(2) * np.array([
            [1],
            [1]
        ])
        
        oop_mat = np.eye(2)
        
        Proj = block_diag(CC_mat,CX_mat,CH_mat,CX_ang_mat,CH_ang_mat,tor_mat,oop_mat)
        
        
        self.Proj = Proj                     
