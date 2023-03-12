import numpy as np
from scipy import linalg
import math
import random

np.set_printoptions(precision=2, suppress=True)


def MBlock(eigenval):
    return np.array(
        (
            [np.cos(eigenval), -np.sin(eigenval)],
            [np.sin(eigenval), np.cos(eigenval)],
        ), dtype=complex
    )


def rot_matrix(matrix):
    # schur decomposition of input matrix
    S, W = linalg.schur(matrix, output='real')

    # block_diag = W.transpose() @ matrix @ W  # block diagonal matrix of eigenvalues
    e = 1j*linalg.eigvals(S)
    m_elements = MBlock(e[0])
    trig_block = np.zeros((len(e), len(e)), dtype=complex)
    # for i in range(2, len(e), 2):
    #     m_elements = np.hstack((m_elements, MBlock(e[i])))
    for i in range(0, len(e), 2):
        trig_block[i, i] = MBlock(e[i])[0, 0]
        trig_block[i, i+1] = MBlock(e[i])[0, 1]
        trig_block[i+1, i] = MBlock(e[i])[1, 0]
        trig_block[i+1, i+1] = MBlock(e[i])[1, 1]
    return W @ trig_block @ W.transpose()


def correlation(L, q_f):
    corr = np.zeros((2*L, 2*L))
    eps = 1 * 10**-10
    for m in range(len(corr)):
        for n in range(len(corr)):
            corr[m, n] = math.sin(q_f*(m - n)) / (math.pi * (m-n + eps))

    return corr


def entang_entropy(correlation_matrix, L):
    """
    Function takes a correlation matrix, bipartitions the matrix, to form M_A. Block diagonalises
    M_A, via a Schur Decomposition. Then finds the eigenvalues of the 2x2 blocks within M_A. Then to
    relate them to the eigenvalues of the reduced density matrix, we use p_i = (\mu_i + 1 )/2, where
    \mu_i are the eigenvalues of the 2x2 blocks. Function returns the Shannon entropy from this, via
    a simple formula. 
    """
    block_size = round(L/2)
    submatrix = correlation_matrix[:2*block_size, :2*block_size]
    M_R, R = linalg.schur(submatrix, output='real')

    for i in range(0, L, 2):
        for j in range(0, L, 2):
            print(linalg.eigvals(M_R[i:i+1, j:j+1]))

    # return M_R, submatrix, linalg.eigvals(M_R)


def H1(site, L):
    H1 = np.zeros((2*L, 2*L))
    H1[site, site+1] = math.pi / 8
    H1[site+1, site] = - math.pi / 8
    return H1

def random_circuit(L):
    R = np.zeros((2*L, 2*L))
    sites = list(range(2*L))
    q = random.sample(sites, 1)[0]
    gate = random.sample(list(range(2)), 1)[0]
    if gate == 1:
        #apply H1
        




# A = np.array(([0, 2, 0, 9], [-2, 0, 6, 4], [0, -6, 0, -1], [-9, -4, 1, 0]))
N = 10
# e = 1j*linalg.eigvals(A)

test_m = correlation(N, math.pi/2)

# print(test_m)


# print(test_m + test_m.transpose())
# print(entang_entropy(test_m, 3)[1])
# print(entang_entropy(test_m, 3)[0])


gate1 = H1(0, 3)
gate2 = H1(2, 3)

print(rot_matrix(gate1))
print(rot_matrix(gate2))
M = rot_matrix(gate1) @ rot_matrix(gate1).transpose()
print(M)

print(entang_entropy(M, 3))
