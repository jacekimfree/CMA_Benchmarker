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

        a, b = 1/4*(-1+5**(1/2)), 1/4*(-1-5**(1/2))
        c, d = (5/8+5**(1/2)/8)**(1/2), (5/8-5**(1/2)/8)**(1/2)

        cyc5_stroop = normalize(np.array([
            [1, 1, 1, 1, 1],
            [1, a, b, b, a],
            [1, b, a, a, b],
            [0, c, d,-d,-c],
            [0, d,-c, c,-d]
        ]).T)

        cyc5_angtor = normalize(np.array([
            [1, b, a, a, b],
            [0, d,-c, c,-d]
        ]).T)

        cyc5_ch_ang = normalize(np.array([
            [1,-1, 1,-1, 1,-1, 1,-1, 1,-1],
            [1,-1, a,-a, b,-b, b,-b, a,-a],
            [1,-1, b,-b, a,-a, a,-a, b,-b],
            [0, 0, c,-c, d,-d,-d, d,-c, c],
            [0, 0, d,-d,-c, c, c,-c,-d, d]
        ]).T)

        Proj = block_diag(cyc5_stroop, cyc5_stroop, cyc5_angtor, cyc5_ch_ang, cyc5_angtor, cyc5_stroop)

        self.Proj = Proj

        # self.sym_sort = np.array([
        #     [],
        #     [],
        #     [],
        #     [],
        # ], dtype=object)

def normalize(mat):
    return 1/norm(mat, axis=0)*mat


if __name__ == "__main__":
    np.set_printoptions(linewidth=400, precision=2, threshold=100000)
    p = Projection([])
    p.run()
    print(p.Proj)
