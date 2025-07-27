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

        a, b = np.cos(45*np.pi/180), np.sin(45*np.pi/180)

        cyc_8str = normalize(np.array([
            [1, 1, 1, 1, 1, 1, 1, 1],
            [1,-1, 1,-1, 1,-1, 1,-1],
            [1, a, 0,-a,-1,-a, 0, a],
            [0, b, 1, b, 0,-b,-1,-b],
            [1, 0,-1, 0, 1, 0,-1, 0],
            [0, 1, 0,-1, 0, 1, 0,-1],
            [1,-a, 0, a,-1, a, 0,-a],
            [0,-b, 1,-b, 0, b,-1, b]
        ]).T)

        ch_8str = normalize(np.array([
            [1, 1, 1, 1, 1, 1, 1, 1],
            [1,-1, 1,-1, 1,-1, 1,-1],
            [1, a, 0,-a,-1,-a, 0, a],
            [0, b, 1, b, 0,-b,-1,-b],
            [1, 0,-1, 0, 1, 0,-1, 0],
            [0, 1, 0,-1, 0, 1, 0,-1],
            [1,-a, 0, a,-1, a, 0,-a],
            [0,-b, 1,-b, 0, b,-1, b]
        ]).T)

        cyc_8ang = normalize(np.array([
            [1,-1, 1,-1, 1,-1, 1,-1],
            [1, 0,-1, 0, 1, 0,-1, 0],
            [0, 1, 0,-1, 0, 1, 0,-1],
            [1,-a, 0, a,-1, a, 0,-a],
            [0,-b, 1,-b, 0, b,-1, b]
        ]).T)

        ch_8ang = normalize(np.array([
            [1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1],
            [1,-1,-1, 1, 1,-1,-1, 1, 1,-1,-1, 1, 1,-1,-1, 1],
            [1,-1, a,-a, 0, 0,-a, a,-1, 1,-a, a, 0, 0, a,-a],
            [0, 0, b,-b, 1,-1, b,-b, 0, 0,-b, b,-1, 1,-b, b],
            [1,-1, 0, 0,-1, 1, 0, 0, 1,-1, 0, 0,-1, 1, 0, 0],
            [0, 0, 1,-1, 0, 0,-1, 1, 0, 0, 1,-1, 0, 0,-1, 1],
            [1,-1,-a, a, 0, 0, a,-a,-1, 1, a,-a, 0, 0,-a, a],
            [0, 0,-b, b, 1,-1,-b, b, 0, 0, b,-b,-1, 1, b,-b]
        ]).T)

        cyc_8tor = normalize(np.array([
            [1,-1, 1,-1, 1,-1, 1,-1],
            [1, 0,-1, 0, 1, 0,-1, 0],
            [0, 1, 0,-1, 0, 1, 0,-1],
            [1,-a, 0, a,-1, a, 0,-a],
            [0,-b, 1,-b, 0, b,-1, b]
        ]).T)

        ch_8oop = normalize(np.array([
            [1, 1, 1, 1, 1, 1, 1, 1],
            [1,-1, 1,-1, 1,-1, 1,-1],
            [1, a, 0,-a,-1,-a, 0, a],
            [0, b, 1, b, 0,-b,-1,-b],
            [1, 0,-1, 0, 1, 0,-1, 0],
            [0, 1, 0,-1, 0, 1, 0,-1],
            [1,-a, 0, a,-1, a, 0,-a],
            [0,-b, 1,-b, 0, b,-1, b]
        ]).T)

        Proj = block_diag(cyc_8str,ch_8str,ch_8ang,cyc_8ang,cyc_8tor,ch_8oop)

        self.Proj = Proj

def normalize(mat):
    return 1/norm(mat, axis=0)*mat


if __name__ == "__main__":
    np.set_printoptions(linewidth=400, precision=2, threshold=100000)
    p = Projection([])
    p.run()
    print(p.Proj)
