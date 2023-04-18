import numpy as np
from scipy import linalg, sparse
import math
import random
import matplotlib.pyplot as plt
np.set_printoptions(precision=3, suppress=True)


def correlation_mat(L):
    # correlation matrix
    # return sparse.diags([1j, 0, -1j], [-1, 0, 1], shape=(2*L, 2*L)).toarray()
    I = np.eye(L)
    pauli_y = np.array(([0, -1j], [1j, 0]))
    return np.kron(I, pauli_y)

# simple gate


def H1(site, L):
    H1 = np.zeros((2*L, 2*L), dtype=complex)
    H1[site, site+1] = math.pi / 8
    H1[site+1, site] = - math.pi / 8
    return H1


def H2(site, L):
    H2 = np.zeros((2*L, 2*L))
    smallh = np.zeros((4, 4))
    smallh[0, 0] = math.pi/8
    smallh[3, 3] = -math.pi/8
    smallh[0+1, 0+2] = math.pi/8
    smallh[0+2, 0+1] = -math.pi / 8
    H2[site:site+4, site:site+4] = smallh
    return H2


def U2(site, L):
    H = np.zeros((2*L, 2*L), dtype=complex)

    smallh = np.zeros((4, 4), dtype=complex)
    smallh[0, 3] = 1
    smallh[3, 0] = -1
    smallh[1, 2] = -1
    smallh[2, 1] = 1
    smallh[0, 1] = 1
    smallh[1, 0] = -1
    smallh[2, 3] = 1
    smallh[3, 2] = -1

    H[site:site+4, site:site+4] = 1*(math.pi)*smallh
    return H


def Xf(site, L):
    X = np.zeros((2*L, 2*L))
    smallh = np.zeros((4, 4))
    smallh[0, 0] = 1
    smallh[3, 3] = -1
    smallh[0+1, 0+2] = -1
    smallh[0+2, 0+1] = 1
    X[site:site+4, site:site+4] = smallh
    return X


def Had(site, L):
    H = np.zeros((2*L, 2*L))
    smallh = np.zeros((4, 4))
    smallh[0, 0] = 1
    smallh[3, 3] = -1
    smallh[1, 2] = 1/np.sqrt(2)
    smallh[2, 1] = 1/np.sqrt(2)
    smallh[1, 1] = 1/np.sqrt(2)
    smallh[2, 2] = -1/np.sqrt(2)
    H[site:site+4, site:site+4] = smallh
    return H


def MBlock(eigenval):
    return np.array(
        (
            [np.cos(eigenval), -np.sin(eigenval)],
            [np.sin(eigenval), np.cos(eigenval)],
        ), dtype=complex
    )


def hop(site, L):
    hop = np.zeros((2*L, 2*L), dtype=complex)
    smallh = np.zeros((4, 4), dtype=complex)
    smallh[0, 2] = 1
    smallh[0, 3] = 1
    smallh[1, 2] = -1
    smallh[1, 3] = 1

    smallh[2, 0] = -1
    smallh[3, 0] = -1
    smallh[2, 1] = 1
    smallh[3, 1] = -1

    # smallh[0, 1] = 1
    # smallh[1, 0] = -1

    # smallh[2, 3] = 1
    # smallh[3, 2] = -1

    hop[site:site+4, site:site+4] = 2*smallh * math.pi
    return hop


def rot_matrix(matrix):
    # schur decomposition of input matrix
    S, W = linalg.schur(matrix, output='real')
    # print(S)
    # print(S)
    e = 1j*linalg.eigvals(S)

    trig_block = np.zeros((len(e), len(e)), dtype=complex)
    for i in range(0, len(e), 2):
        trig_block[i:i+2, i:i+2] = MBlock(e[i])
    # return W @ trig_block @ W.transpose()
    return W@trig_block@W.T


def entang_entropy(correlation_matrix, L):
    """
    Function takes a correlation matrix, bipartitions the matrix, to form M_A. Block diagonalises
    M_A, via a Schur Decomposition. Then finds the eigenvalues of the 2x2 blocks within M_A. Then to
    relate them to the eigenvalues of the reduced density matrix, we use p_i = (\mu_i + 1 )/2, where
    \mu_i are the eigenvalues of the 2x2 blocks. Function returns the Shannon entropy from this, via
    a simple formula.
    """
    Blocksize = round(L/2)

    submatrix = correlation_matrix[:2*Blocksize, 0:2*Blocksize]
    S, R = linalg.schur(correlation_matrix, output='real')
    eigs = []
    # print(S)
    for i in range(0, L, 2):
        eigs.extend(linalg.eigvals(
            submatrix[i:i+2, j:j+2])[0]for j in range(0, L, 2))
    # eigs = linalg.eigvals(submatrix)

    mu = list(filter(lambda num: num > 1*10**(-5), eigs))
    p_list = list(map((lambda mu: (abs(mu) + 1) / 2), mu))
    p = np.asarray(p_list)
    entropy = []
    for i in range(0, len(p), 2):
        if 1-p[i] <= 0:
            entropy.append(-p[i] * math.log(p[i], 2))
        else:
            entropy.append(-p[i] * math.log(p[i], 2) -
                           ((1-p[i]) * math.log((1-p[i]), 2)))

    return sum(entropy)


def random_circuit(L, timesteps):
    sites = list(range(0, L, 2))
    S_A = []
    M = correlation_mat(L)
    for i in range(timesteps):
        q = random.sample(sites, 1)[0]
        gate = random.sample(list(range(2)), 1)[0]
        if gate == 0:
            small_r = rot_matrix(hop(q, L))
        elif gate == 1:
            small_r = rot_matrix(hop(q, L))
        M = np.matmul(small_r, np.matmul(M, small_r.T))

        S_A.append(abs(entang_entropy(M, L)))

    return S_A


####################################################


N = 20
tsteps = 1000
# testgate = H1(0, N)
# test_rotate = rot_matrix(testgate)

test_circ = random_circuit(N, tsteps)
# print(test_circ)
vacuum = correlation_mat(N)
vac_entropy = entang_entropy(vacuum, N)
print(vac_entropy)
plt.plot(range(tsteps), test_circ)
plt.xlabel(r"  t", rotation=0, loc='center')
plt.ylabel(r" $S_{A}/S_{\infty}$      ", rotation=0, loc='center')
# plt.savefig('Fermionic Entanglement Entropy_extended_75.pdf')
plt.show()


# evolve = test_rotate @ vacuum @ test_rotate.T
# print(evolve)
# test_entropy = entang_entropy(evolve, N)
# print(test_entropy)
# def corrtest():
#     M = np.zeros((2*N, 2*N))
#     for i in range(2*N):
#         for j in range(2*N):
#             M[i, j] =


# simple_circ = test_rotate @ vacuum @ test_rotate.T
# print(vac_entropy)


# # entropy = entang_entropy(simple_circ, N)


# # print(test_rotate)
# # print(vac_entropy)
