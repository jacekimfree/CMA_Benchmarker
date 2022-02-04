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

        HA_stretch = np.array([
        1/np.sqrt(3)*np.array([ 1, 1, 1]),
        1/np.sqrt(2)*np.array([ 1,-1, 0]),
        1/np.sqrt(6)*np.array([-1,-1, 2]),
        ]).T
        
        stretch = np.array([
        1/np.sqrt(4)*np.array([ 1, 1, 1, 1]),
        1/np.sqrt(4)*np.array([ 1, 1,-1,-1]),
        1/np.sqrt(4)*np.array([ 1,-1, 1,-1]),
        1/np.sqrt(4)*np.array([ 1,-1,-1, 1])
        ]).T

        angles = np.array([
        0.5/np.sqrt(10)*np.array([ 4, 4,-1,-1,-1,-1,-1,-1,-1,-1]),
        0.5/np.sqrt(10)*np.array([ 4,-4,-1,-1,-1,-1, 1, 1, 1, 1]),
          1/np.sqrt(8)* np.array([ 0, 0, 1, 1,-1,-1, 1, 1,-1,-1]),
          1/np.sqrt(8)* np.array([ 0, 0, 1, 1,-1,-1,-1,-1, 1, 1]),
          1/np.sqrt(8)* np.array([ 0, 0, 1,-1, 1,-1, 1,-1, 1,-1]),
          1/np.sqrt(8)* np.array([ 0, 0, 1,-1, 1,-1,-1, 1,-1, 1]),
          1/np.sqrt(8)* np.array([ 0, 0, 1, 1,-1,-1, 1, 1,-1,-1]),
          1/np.sqrt(8)* np.array([ 0, 0, 1, 1,-1,-1,-1,-1, 1, 1]),
        ]).T

        Proj = block_diag(HA_stretch,stretch,angles,angles2,oop)


        self.Proj = Proj

