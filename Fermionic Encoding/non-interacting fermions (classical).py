import numpy as np
from scipy import linalg, sparse
import math
import random

np.set_printoptions(precision=3, suppress=True)


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
    # print(S)

    # block_diag = W.transpose() @ matrix @ W  # block diagonal matrix of eigenvalues
    e = linalg.eigvals(S)

    # print(e)
    m_elements = MBlock(e[0])
    trig_block = np.zeros((len(e), len(e)))
    # for i in range(2, len(e), 2):
    #     m_elements = np.hstack((m_elements, MBlock(e[i])))
    for i in range(0, len(e), 2):
        # trig_block[i, i] = MBlock(e[i])[0, 0]
        # trig_block[i, i+1] = MBlock(e[i])[0, 1]
        # trig_block[i+1, i] = MBlock(e[i])[1, 0]
        # trig_block[i+1, i+1] = MBlock(e[i])[1, 1]
        trig_block[i:i+2, i:i+2] = MBlock(e[i])
    return W @ trig_block @ W.transpose()


# A = np.array(([0, 2, 0, 9], [-2, 0, 6, 4], [0, -6, 0, -1], [-9, -4, 1, 0]))

# test_rotation = rot_matrix(A)
# print(test_rotation)

#################################################################

# def hopping_correlation(L, q_f):
#     corr = np.zeros((2*L, 2*L))
#     eps = 1 * 10**-10
#     for m in range(len(corr)):
#         for n in range(len(corr)):
#             corr[m, n] = math.sin(q_f*(m - n)) / (math.pi * (m-n + eps))

#     return corr
###########################################################################################


def H1(site, L):
    H1 = np.zeros((2*L, 2*L))
    H1[site, site+1] = math.pi / 8
    H1[site+1, site] = - math.pi / 8
    return H1


# test_gate = H1(0, 8)
# print(test_gate + test_gate.T)


def H2(site, L):
    H2 = np.zeros((2*L, 2*L))
    H2[site, site+3] = math.pi / 2
    H2[3, site] = -math.pi / 2
    H2[site+1, site+2] = -math.pi/2
    H2[site+2, site+1] = math.pi / 2

    return H2


# print(rot_matrix(H2(0, 4)) @ rot_matrix(H2(0, 4)).T)


# # def Xf(site, L):
# #     X = np.zeros((2*L, 2*L))
# #     X[site, site] = 1
# #     X[site+3, site+3] = -1
# #     X[site+1, site+2] = 1
# #     X[site+2, site+1] = 1
# #     return X


# testx = Xf(0, 4)
# print(testx + testx.T)


def logicalHad(site, L):
    Had = np.zeros((2*L, 2*L))
    Had[site, site] = 1
    Had[site+3, site+3] = -1
    Had[site+1, site+1] = 1/np.sqrt(2)
    Had[site+2, site+1] = 1/np.sqrt(2)
    Had[site+1, site+2] = 1/np.sqrt(2)
    Had[site+2, site+2] = -1/np.sqrt(2)
    return Had


def entang_entropy(correlation_matrix, L):
    """
    Function takes a correlation matrix, bipartitions the matrix, to form M_A. Block diagonalises
    M_A, via a Schur Decomposition. Then finds the eigenvalues of the 2x2 blocks within M_A. Then to
    relate them to the eigenvalues of the reduced density matrix, we use p_i = (\mu_i + 1 )/2, where
    \mu_i are the eigenvalues of the 2x2 blocks. Function returns the Shannon entropy from this, via
    a simple formula.
    """

    submatrix = correlation_matrix[:L, :L]
    # print(submatrix)
    M_R, R = linalg.schur(submatrix, output='real')
    # print(M_R)
    eigs = []
    for i in range(0, L, 2):
        eigs.extend(linalg.eigvals(M_R[i:i+1, j:j+1])[0]
                    for j in range(0, L, 2))
    mu = list(filter(lambda num: num != 0, eigs))
    p = list(map((lambda mu: (abs(mu) + 1) / 2), mu))
    good = list(filter(lambda num: 1 - 10e-16 < num < 1 + 10e-16, p))
    entropy = []
    for num in good:
        entropy.append(-num * math.log(num, 2) - (1-num)*math.log((1-num), 2))
        # if num > 1 - 1e-15:
        #     entropy.append(-num * math.log(num, 2))
        # else:
    return sum(entropy)


def correlation_mat(L):
    # return sparse.diags([1j, 1, -1j], [-1, 0, 1], shape=(2*L, 2*L)).toarray()
    I = np.eye(L)
    pauli_y = np.array(([0, -1j], [1j, 0]))
    return np.kron(I, pauli_y)


def random_circuit(L, timesteps):
    R = np.eye(2*L)
    sites = list(range(L))

    for _ in range(timesteps):
        q = random.sample(sites, 1)[0]
        gate = random.sample(list(range(3)), 1)[0]
        if gate == 0:
            small_r = rot_matrix(H1(q, L))
        elif gate == 1:
            small_r = rot_matrix(Xf(q, L))
        elif gate == 2:
            small_r = rot_matrix(logicalHad(q, L))
        R = R@small_r
    return R


# print(rot_matrix(logicalHad(0, 2)))


N = 10
tsteps = 1000
# H1 = logicalHad(3, N)
# H2 = logicalHad(5, N)
# test_circuit = correlation_mat(N)

test = random_circuit(N, tsteps)
print(abs(test).max())
normed = test


# def testy():
#     M = np.empty((2*N, 2*N), dtype=complex)
#     for i in range(2*N):
#         for j in range(2*N):
#             M[i, j] = normed[i, j]*1j*correlation_mat(N)[i, j] * normed[i, j].T

#     return M


# print(testy())
randomcirc = normed @ correlation_mat(N) @ normed.T
# @ correlation_mat(N) @ random_circuit(N, tsteps).transpose()
# print(randomcirc)
entropy_vac = entang_entropy(correlation_mat(N), N)
entrop = entang_entropy(randomcirc, N)
print("Vacuum state entropy: ", entropy_vac)
print("kicked system:", entrop)
