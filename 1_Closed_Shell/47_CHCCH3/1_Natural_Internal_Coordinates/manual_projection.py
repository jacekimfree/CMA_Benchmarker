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

        stretches_mat = np.array([
            [1,  0, 0, 0           ,  0           ,  0           ],
            [0,  1, 0, 0           ,  0           ,  0           ],
            [0,  0, 1, 0           ,  0           ,  0           ],        
            [0,  0, 0, 1/np.sqrt(3),  0           ,  2/np.sqrt(6)],
            [0,  0, 0, 1/np.sqrt(3),  1/np.sqrt(2), -1/np.sqrt(6)],
            [0,  0, 0, 1/np.sqrt(3), -1/np.sqrt(2), -1/np.sqrt(6)]
        ])
        
        ang_mat = np.array([
            [1/np.sqrt(6),  2/np.sqrt(6),  0           ,  0           ,  0           ],
            [1/np.sqrt(6), -1/np.sqrt(6),  1/np.sqrt(2),  0           ,  0           ],
            [1/np.sqrt(6), -1/np.sqrt(6), -1/np.sqrt(2),  0           ,  0           ],
            [-1/np.sqrt(6), 0           ,  0           ,  2/np.sqrt(6),  0           ],
            [-1/np.sqrt(6), 0           ,  0           , -1/np.sqrt(6),  1/np.sqrt(2)],
            [-1/np.sqrt(6), 0           ,  0           , -1/np.sqrt(6), -1/np.sqrt(2)]
        ])
        
        lin_mat = np.eye(4)
        
        Proj = block_diag(stretches_mat,ang_mat,lin_mat)
        
        
        self.Proj = Proj                     
