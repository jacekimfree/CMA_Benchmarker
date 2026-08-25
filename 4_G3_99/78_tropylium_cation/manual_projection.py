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

        a, b, c = np.sin(np.pi/14), np.sin(np.pi/7), np.sin(3*np.pi/14)
        d, e, f = np.cos(np.pi/14), np.cos(np.pi/7), np.cos(3*np.pi/14)

        # 0-6, 7-13, 
        cyc7_str = normalize(np.array([
            [1, 1, 1, 1, 1, 1, 1],  # A1'
            [1, c,-a,-e,-e,-a, c],  # E1'
            [0, f, d, b,-b,-d,-f],
            [1,-a,-e, c, c,-e,-a],  # E2'
            [0, d,-b,-f, f, b,-d],
            [1,-e, c,-a,-a, c,-e],  # E3'
            [0, b,-f, d,-d, f,-b]
        ]).T)

        # 14-17
        cyc7_ang = normalize(np.array([
            [1,-a,-e, c, c,-e,-a],  # E2'
            [0, d,-b,-f, f, b,-d],
            [1,-e, c,-a,-a, c,-e],  # E3'
            [0, b,-f, d,-d, f,-b]
        ]).T)

        # 18-24
        cyc7_ch_ang = normalize(np.array([
            [1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1], # A1'
            [1,-1, c,-c,-a, a,-e, e,-e, e,-a, a, c,-c], # E1'
            [0, 0, f,-f, d,-d, b,-b,-b, b,-d, d,-f, f],
            [1,-1,-a, a,-e, e, c,-c, c,-c,-e, e,-a, a], # E2'
            [0, 0, d,-d,-b, b,-f, f, f,-f, b,-b,-d, d],
            [1,-1,-e, e, c,-c,-a, a,-a, a, c,-c,-e, e], # E3'
            [0, 0, b,-b,-f, f, d,-d,-d, d, f,-f,-b, b]
        ]).T)

        # 25-28
        cyc7_tor = normalize(np.array([
            [1,-a,-e, c, c,-e,-a],  # E2''
            [0, d,-b,-f, f, b,-d],
            [1,-e, c,-a,-a, c,-e],  # E3''
            [0, b,-f, d,-d, f,-b]
        ]).T)

        # 29-35
        cyc7_oop = normalize(np.array([
            [1, 1, 1, 1, 1, 1, 1],  # A1''
            [1, c,-a,-e,-e,-a, c],  # E1''
            [0, f, d, b,-b,-d,-f],
            [1,-a,-e, c, c,-e,-a],  # E2''
            [0, d,-b,-f, f, b,-d],
            [1,-e, c,-a,-a, c,-e],  # E3''
            [0, b,-f, d,-d, f,-b]
        ]).T)

        Proj = block_diag(cyc7_str, cyc7_str, cyc7_ang, cyc7_ch_ang, cyc7_tor, cyc7_oop)

        self.Proj = Proj

        self.sym_sort = np.array([
            [0, 7, 18],
            [],
            [1, 2, 8, 9, 19, 20],
            [3, 4, 10, 11, 14, 15, 21, 22],
            [5, 6, 12, 13, 16, 17, 23, 24],
            [29],
            [],
            [30, 31],
            [25, 26, 32, 33],
            [27, 28, 34, 35]
        ], dtype=object)

def normalize(mat):
    return 1/norm(mat, axis=0)*mat


if __name__ == "__main__":
    np.set_printoptions(linewidth=400, precision=2, threshold=100000)
    p = Projection([])
    p.run()
    print(p.Proj)
