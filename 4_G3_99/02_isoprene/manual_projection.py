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

        ch_2str = normalize(np.array([
            [1, 1],
            [1, -1]
        ]).T)

        ch3_str = normalize(np.array([
            [1, 1, 1],
            [2, -1, -1],
            [0, 1, -1]
        ]).T)

        cc3_plain_ang = normalize(np.array([
            [2, -1, -1],
            [0, 1, -1]
        ]).T)

        ch_ang = normalize(np.array([
            [1, -1]
        ]).T)

        ch2c_ang = normalize(np.array([
            [2, -1, -1],
            [0, 1, -1]
        ]).T)

        ch3_ang = normalize(np.array([
            [1, 1, 1, -1, -1, -1],
            [2, -1, -1, 0, 0, 0],
            [0, 1, -1, 0, 0, 0],
            [0, 0, 0, 2, -1, -1],
            [0, 0, 0, 0, 1, -1]
        ]).T)

        ccc3_rot = normalize(np.array([
            [1, 1]
        ]).T)

        ch2_rot = normalize(np.array([
            [1, 1]
        ]).T)

        ch3_rot = normalize(np.array([
            [1, 1, 1]
        ]).T)

        cc3_plain_oop = normalize(np.array([
            [1, 1, 1]
        ]).T)

        Proj = block_diag(ch3_str, unc, unc, ch_2str, ch_2str, ch3_str, unc, cc3_plain_ang, ch_ang,
                          ch2c_ang, ch2c_ang, ch3_ang, ccc3_rot, ch2_rot, ch2_rot, ch3_rot, cc3_plain_oop, unc, unc, unc)

        projBuff = Proj.copy()
        projBuff = np.append(projBuff[:,:14],oop,axis=1)
        Proj = np.append(projBuff,Proj[:,14:-1],axis=1)
        self.sym_sort = np.array([
            [0, 1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 13,],
            [4, 12, 14, 15, 16, 17, 18]
            ],dtype=object)

        self.Proj = Proj
        np.set_printoptions(edgeitems=60,linewidth=1000)

def normalize(mat):
    return 1/norm(mat, axis=0)*mat


if __name__ == "__main__":
    np.set_printoptions(linewidth=400, precision=2, threshold=100000)
    p = Projection([])
    p.run()
    print(p.Proj)
