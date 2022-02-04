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

        stretch = np.eye(5)

        angles = np.array([
        1/np.sqrt(6)*np.array([ 2,-1,-1]),
        1/np.sqrt(2)*np.array([ 0, 1,-1]),
        ]).T

        tor = np.eye(2)

        #raise RuntimeError

        Proj = block_diag(stretch,angles,tor)


        self.Proj = Proj

