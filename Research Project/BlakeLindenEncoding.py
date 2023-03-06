from scipy import linalg, sparse
import numpy as np
import matplotlib.pyplot as ply
import random
from numpy.linalg import matrix_rank


class Check:

    def __init__(self, N):
        self.N = N  # number of qubits in spin chain

    # def gen_mat(self):
    #     x_string = np.fliplr(sparse.diags(
    #         [0, 1, 0], [-1, 0, 1], shape=(self.N, self.N)).toarray())
    #     return np.concatenate((x_string, np.zeros((self.N, self.N))), axis=1)
    def gen_xmat(self):
        return np.zeros((self.N, self.N))

    def gen_zmat(self):
        return sparse.diags([0, 1, 0], [-1, 0, 1], shape=(self.N, self.N)).toarray()


def tgate(Xmat, Zmat, j):
    temp = Xmat[j].copy()
    Xmat[j] = Zmat[j]
    Zmat[j] = temp
    return Xmat, Zmat


def c_three(Xmat, Zmat, q):
    Zmat[q] = (Zmat[q] + Xmat[q+1] + Zmat[q+1] + Xmat[q+2] + Zmat[q+2]) % 2
    Xmat[q+1] = (Xmat[q] + Xmat[q+1]) % 2
    Zmat[q+1] = (Zmat[q+1] + Xmat[q]) % 2
    Xmat[q+2] = (Xmat[q] + Xmat[q+2]) % 2
    Zmat[q+2] = (Zmat[q+2] + Xmat[q]) % 2
    return Xmat, Zmat


def entangle(N, timesteps):
    BAL = Check(N+1)
    lx = BAL.gen_xmat()
    rz = BAL.gen_zmat()
    S_A = np.zeros(timesteps)
    p = 30

    for i in range(timesteps):
        qubits = list(range(1, N + 1))
        q = random.sample(qubits, 1)[0]
        lx, lz = tgate(lx, rz, q)
        del qubits[N - 3:]
        k = random.sample(qubits, 1)[0]
        lx, lz = c_three(lx, lz, k)
        S_A[i] = matrix_rank(lx[:2*p], 1)

    return lx, lz, S_A


print(entangle(120, 5000)[2])
