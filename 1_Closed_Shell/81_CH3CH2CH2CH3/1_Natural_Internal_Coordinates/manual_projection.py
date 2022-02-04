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
        1/np.sqrt(2)*np.array([ 1, 1, 0]),
        1/np.sqrt(2)*np.array([ 1,-1, 0]),
        1/np.sqrt(1)*np.array([ 0, 0, 1]),
        ]).T
        
        inner_CH_str = np.array([
        1/np.sqrt(4)*np.array([ 1, 1, 1, 1]),
        1/np.sqrt(4)*np.array([ 1, 1,-1,-1]),
        1/np.sqrt(4)*np.array([ 1,-1, 1,-1]),
        1/np.sqrt(4)*np.array([ 1,-1,-1, 1]),
        ]).T

        outer_CH_str = np.array([
         1/np.sqrt(6)*np.array([1, 1, 1, 1, 1, 1]),
         1/np.sqrt(6)*np.array([1, 1, 1,-1,-1,-1]),
        1/np.sqrt(12)*np.array([2,-1,-1, 2,-1,-1]),
        1/np.sqrt(12)*np.array([2,-1,-1,-2, 1, 1]),
         1/np.sqrt(2)*np.array([0, 1,-1, 0, 1,-1]),
         1/np.sqrt(2)*np.array([0, 1,-1, 0,-1, 1]),
        ]).T

        HA_ang = np.array([
        1/np.sqrt(2)*np.array([1, 1]),
        1/np.sqrt(2)*np.array([1,-1]),
        ]).T

        inner_ang = np.array([
        0.5/np.sqrt(10)*np.array([4,-1,-1,-1,-1, 4,-1,-1,-1,-1]),
        0.5/np.sqrt(10)*np.array([4,-1,-1,-1,-1,-4, 1, 1, 1, 1]),
          1/np.sqrt(10)*np.array([0, 1, 1,-1,-1, 0, 1, 1,-1,-1]),
          1/np.sqrt(10)*np.array([0, 1, 1,-1,-1, 0,-1,-1, 1, 1]),
          1/np.sqrt(10)*np.array([0, 1,-1, 1,-1, 0, 1,-1, 1,-1]),
          1/np.sqrt(10)*np.array([0, 1,-1, 1,-1, 0,-1, 1,-1, 1]),
          1/np.sqrt(10)*np.array([0, 1,-1,-1, 1, 0, 1,-1,-1, 1]),
          1/np.sqrt(10)*np.array([0, 1,-1,-1, 1, 0,-1, 1, 1,-1]),
        ]).T

        outer_ang = np.array([
        1/np.sqrt(12)*np.array([1, 1, 1,-1,-1,-1, 1, 1, 1,-1,-1,-1]),
        1/np.sqrt(12)*np.array([1, 1, 1,-1,-1,-1,-1,-1,-1, 1, 1, 1]),
         1/np.sqrt(8)*np.array([2,-1,-1, 0, 0, 0, 2,-1,-1, 0, 0, 0]),
         1/np.sqrt(8)*np.array([2,-1,-1, 0, 0, 0,-2, 1, 1, 0, 0, 0]),
         1/np.sqrt(4)*np.array([0, 1,-1, 0, 0, 0, 0, 1,-1, 0, 0, 0]),
         1/np.sqrt(4)*np.array([0, 1,-1, 0, 0, 0, 0,-1, 1, 0, 0, 0]),
         1/np.sqrt(8)*np.array([0, 0, 0, 2,-1,-1, 0, 0, 0, 2,-1,-1]),
         1/np.sqrt(8)*np.array([0, 0, 0, 2,-1,-1, 0, 0, 0,-2, 1, 1]),
         1/np.sqrt(4)*np.array([0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1]),
         1/np.sqrt(4)*np.array([0, 0, 0, 0, 1,-1, 0, 0, 0, 0,-1, 1]),
        ]).T

        tor = np.array([
        1/np.sqrt(1)*np.array([1, 0, 0, 0, 0, 0, 0]),
        1/np.sqrt(6)*np.array([0, 1, 1, 1, 1, 1, 1]),
        1/np.sqrt(6)*np.array([0, 1, 1, 1,-1,-1,-1]),
        ]).T

        Proj = block_diag(HA_stretch,inner_CH_str,outer_CH_str,HA_ang,inner_ang,outer_ang,tor)


        self.Proj = Proj

