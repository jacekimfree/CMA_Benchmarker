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

        # 0-1, 2-3, 4-5, 6-7
        cc_2str = normalize(np.array([
            [1, 1],
            [1,-1]
        ]).T)

        # 8-13
        ch3_2str = normalize(np.array([
            [1, 1, 1, 1, 1, 1],
            [1, 1, 1,-1,-1,-1],
            [2,-1,-1, 2,-1,-1],
            [2,-1,-1,-2, 1, 1],
            [0, 1,-1, 0, 1,-1],
            [0, 1,-1, 0,-1, 1]
        ]).T)

        # 14-15
        ch2_str = normalize(np.array([
            [1, 1],
            [1,-1]
        ]).T)

        # 16-19, 20-23, 24-27
        ch2_2str = normalize(np.array([
            [1, 1, 1, 1],
            [1, 1,-1,-1],
            [1,-1, 1,-1],
            [1,-1,-1, 1]
        ]).T)

        # 28
        unc = np.eye(1)

        # 29-30, 31-32, 33-34
        cc_2ang = normalize(np.array([
            [1, 1],
            [1,-1]
        ]).T)

        # 35-38
        ch2_ang = normalize(np.array([
            [4,-1,-1,-1,-1],
            [0, 1, 1,-1,-1],
            [0, 1,-1, 1,-1],
            [0, 1,-1,-1, 1]
        ]).T)

        # 39-46, 47-54, 55-62
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

        # 63-72
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

        # 73-74, 75-76, 77-78
        cc_2tor = normalize(np.array([
            [1, 1],
            [1,-1]
        ]).T)

        ch3_2rot = normalize(np.array([
            [1, 1, 1, 1, 1, 1],
            [1, 1, 1,-1,-1,-1]
        ]).T)

        Proj = block_diag(cc_2str,cc_2str,cc_2str,cc_2str,ch3_2str,ch2_str,ch2_2str,ch2_2str,ch2_2str,unc,cc_2ang,cc_2ang,cc_2ang,ch2_ang,ch2_2ang,ch2_2ang,ch2_2ang,ch3_2ang,cc_2tor,cc_2tor,cc_2tor,ch3_2rot)

        self.Proj = Proj

        self.sym_sort = np.array([
            [0, 2, 4, 6, 8, 10, 14, 16, 20, 24, 28, 29, 31, 33, 35, 39, 43, 47, 51, 55, 59, 63, 65, 69],
            [13, 19, 23, 27, 38, 42, 46, 50, 54, 58, 62, 68, 72, 73, 75, 77, 79],
            [12, 15, 18, 22, 26, 36, 41, 45, 49, 53, 57, 61, 67, 71, 74, 76, 78, 80],
            [1, 3, 5, 7, 9, 11, 17, 21, 25, 30, 32, 34, 37, 40, 44, 48, 52, 56, 60, 64, 66, 70],
        ], dtype=object)

def normalize(mat):
    return 1/norm(mat, axis=0)*mat


if __name__ == "__main__":
    np.set_printoptions(linewidth=400, precision=2, threshold=100000)
    p = Projection([])
    p.run()
    print(p.Proj)
