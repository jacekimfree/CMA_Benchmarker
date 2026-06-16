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

        # 11, 18, 40, 47
        unc = np.eye(1)

        a, x = np.sin(3*np.pi/14), np.cos(3*np.pi/14)
        b, y = np.sin(np.pi/14), np.cos(np.pi/14)
        c, z = np.sin(np.pi/7), np.cos(np.pi/7)

        # 0-6
        cyc7_str = normalize(np.array([
            [1, 1, 1, 1, 1, 1, 1],
            [1, a,-b,-z,-z,-b, a],
            [0, x, y, c,-c,-y,-x],
            [1,-b,-z, a, a,-z,-b],
            [0, y,-c,-x, x, c,-y],
            [1,-z, a,-b,-b, a,-z],
            [0, c,-x, y,-y, x,-c]
        ]).T)

        # 7-8, 9-10
        cc_2str = normalize(np.array([
            [1, 1],
            [1,-1]
        ]).T)

        # 12-13, 14-15, 16-17
        ch_2str = normalize(np.array([
            [1, 1],
            [1,-1]
        ]).T)

        a, b = np.cos(4*np.pi/5), np.cos(8*np.pi/5)
        c, d = np.sin(4*np.pi/5), np.sin(8*np.pi/5)

        # 19-20, 33-34
        cyc5_angtor = normalize(np.array([
            [1, a, b, b, a],
            [0, c, d,-d,-c]
        ]).T)

        a, b, c = -np.sin(np.pi/14), -np.cos(np.pi/7), np.sin(3*np.pi/14)
        d, e, f = np.cos(np.pi/14), -np.sin(np.pi/7), -np.cos(3*np.pi/14)

        # 21-24, 35-38
        cyc7_angtor = normalize(np.array([
            [1, a, b, c, c, b, a],
            [1, b, c, a, a, c, b],
            [0, d, e, f,-f,-e,-d],
            [0,-e, f, d,-d,-f, e]
        ]).T)

        # 25, 32
        ch_ang = normalize(np.array([
            [1,-1]
        ]).T)

        # 26-27, 28-29, 30-31
        ch_2ang = normalize(np.array([
            [1,-1, 1,-1],
            [1,-1,-1, 1]
        ]).T)

        # 39
        fus2_but = normalize(np.array([
            [1,-1]
        ]).T)

        # 41-42, 43-44, 45-46
        ch_2oop = normalize(np.array([
            [1, 1],
            [1,-1],
        ]).T)

        Proj = block_diag(cyc7_str, cc_2str, cc_2str, unc, ch_2str, ch_2str, ch_2str, unc, cyc5_angtor, cyc7_angtor, ch_ang, ch_2ang, ch_2ang, ch_2ang, ch_ang, cyc5_angtor, cyc7_angtor, fus2_but, unc, ch_2oop, ch_2oop, ch_2oop, unc)

        self.Proj = Proj

        self.sym_sort = np.array([
            [0, 1, 3, 5, 7, 9, 11, 12, 14, 16, 18, 19, 21, 22, 26, 28, 30],
            [33, 35, 36, 41, 43, 45],
            [34, 37, 38, 39, 40, 42, 44, 46, 47],
            [2, 4, 6, 8, 10, 13, 15, 17, 20, 23, 24, 25, 27, 29, 31, 32],
        ], dtype=object)

def normalize(mat):
    return 1/norm(mat, axis=0)*mat


if __name__ == "__main__":
    np.set_printoptions(linewidth=400, precision=2, threshold=100000)
    p = Projection([])
    p.run()
    print(p.Proj)
