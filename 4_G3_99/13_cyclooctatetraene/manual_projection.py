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

        unc = np.eye(1)

        cc_8str = normalize(np.array([
            [1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1,-1,-1,-1,-1],
            [1, 1,-1,-1, 1, 1,-1,-1],
            [1, 1,-1,-1,-1,-1, 1, 1],
            [1,-1, 1,-1, 1,-1, 1,-1],
            [1,-1, 1,-1,-1, 1,-1, 1],
            [1,-1,-1, 1, 1,-1,-1, 1],
            [1,-1,-1, 1,-1, 1, 1,-1],
        ]).T)

        ch_8str = normalize(np.array([
            [1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1,-1,-1,-1,-1],
            [1, 1,-1,-1, 1, 1,-1,-1],
            [1, 1,-1,-1,-1,-1, 1, 1],
            [1,-1, 1,-1, 1,-1, 1,-1],
            [1,-1, 1,-1,-1, 1,-1, 1],
            [1,-1,-1, 1, 1,-1,-1, 1],
            [1,-1,-1, 1,-1, 1, 1,-1],
        ]).T)

        cyc_8ang = np.eye(1)

        ch_8ang = np.eye(1)

        cyc_8tor = np.eye(1)

        ch_8oop = normalize(np.array([
            [1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1,-1,-1,-1,-1],
            [1, 1,-1,-1, 1, 1,-1,-1],
            [1, 1,-1,-1,-1,-1, 1, 1],
            [1,-1, 1,-1, 1,-1, 1,-1],
            [1,-1, 1,-1,-1, 1,-1, 1],
            [1,-1,-1, 1, 1,-1,-1, 1],
            [1,-1,-1, 1,-1, 1, 1,-1],
        ]).T)

        Proj = block_diag(cc_8str,ch_8str,ch_8ang,cyc_8ang,cyc_8tor,ch_8oop)

        self.Proj = Proj
        '''
        self.sym_sort = np.array([                   
    		[],
    		[]              
    		],dtype=object)                          
	'''


def normalize(mat):
    return 1/norm(mat, axis=0)*mat


if __name__ == "__main__":
    np.set_printoptions(linewidth=400, precision=2, threshold=100000)
    p = Projection([])
    p.run()
    print(p.Proj)
