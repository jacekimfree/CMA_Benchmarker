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

        cc_2str = normalize(np.array([
            [1, 1],
            [1,-1]
        ]).T)

        ch2_2str = normalize(np.array([
            [1, 1, 1, 1],
            [1, 1,-1,-1],
            [1,-1, 1,-1],
            [1,-1,-1, 1]
        ]).T)

        ch3_2str = normalize(np.array([
            [1, 1, 1, 1, 1, 1],
            [1, 1, 1,-1,-1,-1],
            [2,-1,-1, 2,-1,-1],
            [2,-1,-1,-2, 1, 1],
            [0, 1,-1, 0, 1,-1],
            [0, 1,-1, 0,-1, 1]
        ]).T)

        cc_2ang = normalize(np.array([
            [1, 1],
            [1,-1]
        ]).T)

        ch2_2ang = normalize(np.array([
            [4,-1,-1,-1,-1, 4,-1,-1,-1,-1],
            [4,-1,-1,-1,-1,-4, 1, 1, 1, 1],
            [0, 1, 1,-1,-1, 0, 1, 1,-1,-1],
            [0, 1, 1,-1,-1, 0,-1,-1, 1, 1],
            [0, 1,-1, 1,-1, 0, 1,-1, 1,-1],
            [0, 1,-1, 1,-1, 0,-1, 1,-1, 1],
            [0, 1,-1,-1, 1, 0, 1,-1,-1, 1],
            [0, 1,-1,-1, 1, 0,-1, 1, 1,-1]
        ]).T)

        ch3_2ang = normalize(np.array([
            [1, 1, 1,-1,-1,-1, 1, 1, 1,-1,-1,-1],
            [1, 1, 1,-1,-1,-1,-1,-1,-1, 1, 1, 1],
            [2,-1,-1, 0, 0, 0, 2,-1,-1, 0, 0, 0],
            [2,-1,-1, 0, 0, 0,-2, 1, 1, 0, 0, 0],
            [0, 1,-1, 0, 0, 0, 0, 1,-1, 0, 0, 0],
            [0, 1,-1, 0, 0, 0, 0,-1, 1, 0, 0, 0],
            [0, 0, 0, 2,-1,-1, 0, 0, 0, 2,-1,-1],
            [0, 0, 0, 2,-1,-1, 0, 0, 0,-2, 1, 1],
            [0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1],
            [0, 0, 0, 0, 1,-1, 0, 0, 0, 0,-1, 1]
        ]).T)

        cc_2tor = normalize(np.array([
            [1, 1],
            [1,-1]
        ]).T)

        ch3_2rot = normalize(np.array([
            [1, 1, 1, 1, 1, 1],
            [1, 1, 1,-1,-1,-1]
        ]).T)

        Proj = block_diag(unc,cc_2str,cc_2str,cc_2str,cc_2str,ch2_2str,ch2_2str,ch2_2str,ch2_2str,ch3_2str,cc_2ang,cc_2ang,cc_2ang,cc_2ang,ch2_2ang,ch2_2ang,ch2_2ang,ch2_2ang,ch3_2ang,unc,cc_2tor,cc_2tor,cc_2tor,ch3_2rot)

        self.Proj = Proj

        self.sym_sort = np.array([
            [0, 1, 3, 5, 7, 9, 13, 17, 21, 25, 27, 31, 33, 35, 37, 39, 43, 47, 51, 55, 59, 63, 67, 71, 73, 77],
            [12, 16, 20, 24, 30, 42, 46, 50, 54, 58, 62, 66, 70, 76, 80, 83, 85, 87, 89],
            [11, 15, 19, 23, 29, 41, 45, 49, 53, 57, 61, 65, 69, 75, 79, 81, 82, 84, 86, 88],
            [2, 4, 6, 8, 10, 14, 18, 22, 26, 28, 32, 34, 36, 38, 40, 44, 48, 52, 56, 60, 64, 68, 72, 74, 78],
        ], dtype=object)

def normalize(mat):
    return 1/norm(mat, axis=0)*mat


if __name__ == "__main__":
    np.set_printoptions(linewidth=400, precision=2, threshold=100000)
    p = Projection([])
    p.run()
    print(p.Proj)
