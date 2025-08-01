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

        cc_6str = normalize(np.array([
            [1, 1, 1, 1, 1, 1],
            [1, -1, 1, -1, 1, -1],
            [2, 1, -1, -2, -1, 1],
            [0, 1, 1, 0, -1, -1],
            [2, -1, -1, 2, -1, -1],
            [0, 1, -1, 0, 1, -1]
        ]).T)

        cc_2str = normalize(np.array([
            [1, 1],
            [1, -1]
        ]).T)

        cc_4str = normalize(np.array([
            [1, 1, 1, 1],
            [1, 1, -1, -1],
            [1, -1, 1, -1],
            [1, -1, -1, 1]
        ]).T)

        cyc_6ang = normalize(np.array([
            [1, -1, 1, -1, 1, -1],
            [2, -1, -1, 2, -1, -1],
            [0, 1, -1, 0, 1, -1]
        ]).T)

        ch_2ang = normalize(np.array([
            [1, -1, 1, -1],
            [1, -1, -1, 1]
        ]).T)

        ch_4ang = normalize(np.array([
            [1, -1, 1, -1, 1, -1, 1, -1],
            [1, -1, 1, -1, -1, 1, -1, 1],
            [1, -1, -1, 1, 1, -1, -1, 1],
            [1, -1, -1, 1, -1, 1, 1, -1]
        ]).T)

        cyc_6tor = normalize(np.array([
            [1, -1, 1, -1, 1, -1],
            [0, -1, 1, 0, -1, 1],
            [2, -1, -1, 2, -1, -1]
        ]).T)

        cc_2tor = normalize(np.array([
            [1, 1],
            [1, -1]
        ]).T)

        ch_2oop = normalize(np.array([
            [1, 1],
            [1, -1],
        ]).T)

        ch_4oop_sym = normalize(np.array([
            [1, 1, 1, 1],
            [1, 1, -1, -1]
        ]).T)

        Proj = block_diag(cc_6str, cc_2str, cc_4str, cyc_6ang, ch_2ang,
                          ch_4ang, cyc_6tor, cc_2tor, ch_2oop, ch_4oop_sym)

        self.Proj = Proj

        self.sym_sort = np.array([
            [0, 4, 6, 8, 13, 17],
            [29],
            [21, 25, 27],
            [5, 11, 14, 16, 20],
            [23, 24],
            [3, 7, 10, 12, 18],
            [1, 2, 9, 15, 19],
            [22, 26, 28],
        ], dtype=object)

def normalize(mat):
    return 1/norm(mat, axis=0)*mat


if __name__ == "__main__":
    np.set_printoptions(linewidth=400, precision=2, threshold=100000)
    p = Projection([])
    p.run()
    print(p.Proj)
