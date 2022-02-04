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

        HC_str = np.eye(1)
 
        HA_str = np.array([
        1/np.sqrt(3)*np.array([ 1, 1, 1]),
        1/np.sqrt(6)*np.array([ 2,-1,-1]),
        1/np.sqrt(2)*np.array([ 0, 1,-1]),
        ]).T
       
        CH_str = np.array([
        1/np.sqrt( 9)*np.array([1, 1, 1, 1, 1, 1, 1, 1, 1]),
        1/np.sqrt(18)*np.array([2, 2, 2,-1,-1,-1,-1,-1,-1]),
        1/np.sqrt( 6)*np.array([0, 0, 0, 1, 1, 1,-1,-1,-1]),
        1/np.sqrt(18)*np.array([2,-1,-1, 2,-1,-1, 2,-1,-1]),
        1/np.sqrt(36)*np.array([4,-2,-2,-2, 1, 1,-2, 1, 1]),
        1/np.sqrt(12)*np.array([0, 0, 0, 2,-1,-1,-2, 1, 1]),
        1/np.sqrt( 6)*np.array([0, 1,-1, 0, 1,-1, 0, 1,-1]),
        1/np.sqrt(12)*np.array([0, 2,-2, 0,-1, 1, 0,-1, 1]),
        1/np.sqrt( 4)*np.array([0, 0, 0, 0, 1,-1, 0,-1, 1]),
        ]).T

        HC_ang = np.array([
        1/np.sqrt(6)*np.array([2,-1,-1]),
        1/np.sqrt(2)*np.array([0, 1,-1]),
        ]).T

        HA_ang = np.array([
        1/np.sqrt(6)*np.array([2,-1,-1]),
        1/np.sqrt(3)*np.array([0, 1,-1]),
        ]).T

        CH_ang = np.array([
        1/np.sqrt(18)*np.array([1, 1, 1,-1,-1,-1, 1, 1, 1,-1,-1,-1, 1, 1, 1,-1,-1,-1]),
        1/np.sqrt(36)*np.array([2, 2, 2,-2,-2,-2,-1,-1,-1, 1, 1, 1,-1,-1,-1, 1, 1, 1]),
        1/np.sqrt(12)*np.array([0, 0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1,-1,-1,-1, 1, 1, 1]),
        1/np.sqrt(18)*np.array([2,-1,-1, 0, 0, 0, 2,-1,-1, 0, 0, 0, 2,-1,-1, 0, 0, 0]),
        1/np.sqrt(36)*np.array([4,-2,-2, 0, 0, 0,-2, 1, 1, 0, 0, 0,-2, 1, 1, 0, 0, 0]),
        1/np.sqrt(12)*np.array([0, 0, 0, 0, 0, 0, 2,-1,-1, 0, 0, 0,-2, 1, 1, 0, 0, 0]),
        1/np.sqrt( 6)*np.array([0, 1,-1, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1, 0, 0, 0]),
        1/np.sqrt(12)*np.array([0, 2,-2, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0,-1, 1, 0, 0, 0]),
        1/np.sqrt( 4)*np.array([0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0,-1, 1, 0, 0, 0]),
        1/np.sqrt(18)*np.array([0, 0, 0, 2,-1,-1, 0, 0, 0, 2,-1,-1, 0, 0, 0, 2,-1,-1]),
        1/np.sqrt(36)*np.array([0, 0, 0, 4,-2,-2, 0, 0, 0,-2, 1, 1, 0, 0, 0,-2, 1, 1]),
        1/np.sqrt(12)*np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-1,-1, 0, 0, 0,-2, 1, 1]),
        1/np.sqrt( 6)*np.array([0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0,-1, 1]),
        1/np.sqrt(12)*np.array([0, 0, 0, 0, 2,-2, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 1,-1]),
        1/np.sqrt( 4)*np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1]),
        ]).T

        tor = np.array([
        1/np.sqrt(18)*np.array([1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1]),
        1/np.sqrt(36)*np.array([2,-2, 2,-2, 2,-2,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1]),
        1/np.sqrt(12)*np.array([0, 0, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1,-1, 1,-1, 1,-1, 1]),
        ]).T 

        oop = np.array([
        1/np.sqrt(3)*np.array([1, 1, 1]),
        ]).T 

        Proj = block_diag(HC_str,HA_str,CH_str,HC_ang,HA_ang,CH_ang,tor,oop)


        self.Proj = Proj

