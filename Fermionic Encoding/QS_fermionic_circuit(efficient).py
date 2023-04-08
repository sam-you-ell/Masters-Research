import numpy as np
import scipy
from scipy.sparse import lil_matrix, csr_matrix, kron, coo_matrix, linalg
import math
import random
import matplotlib.pyplot as plt
import pandas as pd
# np.set_printoptions(precision=2, suppress=True)


def correlation_mat(L):
    # correlation matrix
    # return sparse.diags([1j, 0, -1j], [-1, 0, 1], shape=(2*L, 2*L)).toarray()
    I = lil_matrix(np.eye(L))
    pauli_y = np.array(([0, -1j], [1j, 0]))
    sparse_y = lil_matrix(pauli_y)

    return kron(I, sparse_y)


# simple gate


def H1(site, L):
    H1 = np.zeros((2*L, 2*L), dtype=complex)
    H1[site, site+1] = math.pi/8
    H1[site+1, site] = -math.pi/8
    return H1


def Hodd(site, L):
    H1 = np.zeros((2*L, 2*L), dtype=complex)
    H1[site, site+2] = math.pi/8
    H1[site+2, site] = -math.pi/8
    return H1


def H2(site, L):
    H2 = np.zeros((2*L, 2*L), dtype=complex)
    smallh = np.zeros((4, 4))
    smallh[1, 3] = math.pi/8
    smallh[3, 1] = -math.pi/8
    smallh[2, 3] = math.pi/2
    smallh[3, 2] = -math.pi / 2
    H2[site:site+4, site:site+4] = smallh
    return H2


# def Xf(site, L):
#     X = np.zeros((2*L, 2*L))
#     smallh = np.zeros((4, 4))
#     smallh[0, 0] = 1
#     smallh[3, 3] = -1
#     smallh[0+1, 0+2] = -1
#     smallh[0+2, 0+1] = 1
#     X[site:site+4, site:site+4] = smallh
#     return X


def Had(site, L):
    H = np.zeros((2*L, 2*L), dtype=complex)
    smallh = np.zeros((4, 4))
    smallh[0, 0] = 1
    smallh[3, 3] = -1
    smallh[1, 2] = 1/np.sqrt(2)
    smallh[2, 1] = 1/np.sqrt(2)
    smallh[1, 1] = 1/np.sqrt(2)
    smallh[2, 2] = -1/np.sqrt(2)
    H[site:site+4, site:site+4] = smallh
    return H


# def MBlock(eigenval):
#     return np.array(
#         (
#             [np.cos(eigenval), -np.sin(eigenval)],
#             [np.sin(eigenval), np.cos(eigenval)],
#         ), dtype=complex
#     )

def swap(site, L):
    H = np.zeros((2*L, 2*L), dtype=complex)

    smallh = np.zeros((4, 4))
    smallh[0, 3] = -1
    smallh[3, 0] = 1
    smallh[1, 2] = 1
    smallh[2, 1] = -1
    smallh[0, 1] = 1
    smallh[1, 0] = -1
    smallh[2, 3] = 1
    smallh[3, 2] = -1
    H[site:site+4, site:site+4] = smallh * math.pi/4
    return H


def correctswap(site, L):
    H = np.zeros((2*L, 2*L), dtype=complex)

    smallh = np.zeros((4, 4))
    smallh[0, 0] = 1
    smallh[3, 3] = 1
    smallh[1, 2] = 1
    smallh[2, 1] = -1

    H[site:site+4, site:site+4] = smallh
    return H


# s = swap(0, 4)
# print(s)
# print(H1(0, 2))


def rot_matrix(matrix):
    # schur decomposition of input matrix
    S, W = scipy.linalg.schur(matrix, output='real')
    # print(S)
    # print(S)
    W = lil_matrix(W)

    # evals = scipy.linalg.eigvals(S)
    e = 1j*S.diagonal()

    trig_block = lil_matrix((len(e), len(e)), dtype=complex)
    # trig_block = np.zeros((len(e), len(e)), dtype=complex)
    for i in range(0, len(e), 2):
        trig_block[i, i] = np.cos(e[i])
        trig_block[i, i+1] = -np.sin(e[i])
        trig_block[i+1, i] = np.sin(e[i])
        trig_block[i+1, i+1] = np.cos(e[i])
    return (W@trig_block@W.T)


def entang_entropy(correlation_matrix, L):
    """
    Function takes a correlation matrix, bipartitions the matrix, to form M_A. Block diagonalises
    M_A, via a Schur Decomposition. Then finds the eigenvalues of the 2x2 blocks within M_A. Then to
    relate them to the eigenvalues of the reduced density matrix, we use p_i = (\mu_i + 1 )/2, where
    \mu_i are the eigenvalues of the 2x2 blocks. Function returns the Shannon entropy from this, via
    a simple formula.
    """
    Blocksize = round(L/2)
    corr = correlation_matrix.toarray()
    # Blocksize = round(L/2)
    submatrix = corr[:2*Blocksize, :2*Blocksize]

    # S, R = linalg.schur(submatrix, output='real')
    eigs = []
    # print(S)
    for i in range(0, L, 2):
        eigs.extend(scipy.linalg.eigvals(
            submatrix[i:i+2, j:j+2])[0] for j in range(0, L, 2))
    # eigs = linalg.eigvals(submatrix)

    mu = list(filter(lambda num: num > 1*10**(-5), eigs))
    p_list = list(map((lambda mu: (abs(mu) + 1) / 2), mu))
    p = np.asarray(p_list)
    entropy = []
    for i in range(0, len(p), 2):
        if 1-p[i] <= 0:
            entropy.append(-p[i] * math.log(p[i], 2)/math.sqrt(L))
        else:
            entropy.append((-p[i] * math.log(p[i], 2) -
                            ((1-p[i]) * math.log((1-p[i]), 2)))/math.sqrt(L))

    return sum(entropy)


def random_circuit(L, timesteps):
    sites = list(range(0, L, 1))
    M = correlation_mat(L)
    S_A = [entang_entropy(M, L)]
    for i in range(timesteps-1):
        q = random.sample(sites, 1)[0]
        gate = random.sample(list(range(2)), 1)[0]
        if gate == 0:
            small_r = csr_matrix(rot_matrix(H1(q, L)))
        elif gate == 1:
            small_r = csr_matrix(rot_matrix(H1(q, L)))

        M = (small_r @ M @ small_r.T)
        # for i in range(2*L):
        #     M[i, i] = 0

        S_A.append(abs(entang_entropy(M, L) * math.log(2)/math.log(2**L)))

    return S_A


def commutator(matrix1, matrix2):
    A = matrix1.toarray()
    B = matrix2.toarray()
    return A@B - B@A
####################################################


# N = 50
# tsteps = 4000

N = 30
tsteps = 2000

S = random_circuit(N, tsteps)
t = range(tsteps)

# M = correlation_mat(N)

# print(commutator(S, M))
# Entropy = {
#     "Time": t,
#     "Ent1": random_circuit(N, tsteps),
#     "Ent2": random_circuit(N, tsteps),
#     "Ent3": random_circuit(N, tsteps),
#     "Ent4": random_circuit(N, tsteps),
#     "Ent5": random_circuit(N, tsteps),
# }
# df2 = pd.DataFrame(Entropy)
# with pd.ExcelWriter('FermionCircuit60.xlsx') as writer:
#     df2.to_excel(writer, sheet_name='60_sites')

# print('done')

# print("Page Value:", math.log(2**N))
plt.plot(range(tsteps), S)
plt.xlabel(r"  t", rotation=0, loc='center')
plt.ylabel(r" $S_{A}/S_{\infty}$      ", rotation=0, loc='center')


# # plt.savefig('Fermionic Entanglement Entropy_extended_75.pdf')
plt.show()
