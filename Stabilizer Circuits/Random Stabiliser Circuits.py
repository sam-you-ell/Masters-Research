import numpy as np
from scipy import sparse
import random
import sympy
from numpy.linalg import matrix_rank
import matplotlib.pyplot as plt
import math
import pandas as pd


def checkmatrix(N):
    """
    Function sets up the check matrix, for the simplest case, will likely need to change so I can
    specify a string.
    """
    return sparse.diags([0, 1, 0], [-1, 0, 1], shape=(N, 2*N)).toarray()
    # X = np.zeros((N, N))
    # Z = sparse.diags([0, 1, 0], [-1, 0, 1], shape=(N, N)).toarray()
    # return np.concatenate((X, Z), axis=1)
    # c = np.zeros((N, 2*N))
    # c[3, 3] = 1
    # return c


def Tgate(matrix, q):
    """
    Function to implement the action of the T gate described in the Blake and Linden paper. As the
    paper states, it simply swaps the two rows on either side of a check matrix for a given row
    number and returns the update matrix
    """
    N = len(matrix)
    initialrow = matrix[q, :N].copy()
    matrix[q, :N] = matrix[q, N: 2*N]
    matrix[q, N:2*N] = initialrow
    return matrix


def C3gate(matrix, q, control):
    N = len(matrix)

    if control == q:
        pass
    elif control == q+1:
        matrix[[q, q+1]] = matrix[[q+1, q]]
    elif control == q+2:
        matrix[[q, q+2]] = matrix[[q+2, q]]

    xvec = matrix[:, :N]
    zvec = matrix[:, N: 2*N]
    zvec[q] = (zvec[q] + xvec[q+1] +
               zvec[q+1] + xvec[q+2] + zvec[q+2])
    xvec[q+1] = (xvec[q] + xvec[q+1])
    zvec[q+1] = (zvec[q+1] + xvec[q])
    xvec[q+2] = (xvec[q] + xvec[q+2])
    zvec[q+2] = (zvec[q+2] + xvec[q])

    matrix[:, :N] = xvec % 2
    matrix[:, N: 2*N] = zvec % 2

    if control == q:
        pass
    elif control == q+1:
        matrix[[q+1, q]] = matrix[[q, q+1]]
    elif control == q+2:
        matrix[[q+2, q]] = matrix[[q, q+2]]

    return matrix


def getCutStabilizers(binaryMatrix, cut):
    """
        - Purpose: Return only the part of the binary matrix that corresponds to the qubits we want to consider for a bipartition.
        - Inputs:
            - binaryMatrix (array of size (N, 2N)): The binary matrix for the stabilizer generators.
            - cut (integer): Location for the cut.
        - Outputs:
            - cutMatrix (array of size (N, 2cut)): The binary matrix for the cut on the left.
    """
    N = len(binaryMatrix)
    cutMatrix = np.zeros((N, 2*cut))

    cutMatrix[:, :cut] = binaryMatrix[:, :cut]
    cutMatrix[:, cut:] = binaryMatrix[:, N:N+cut]

    return cutMatrix


def randomcircuit(matrix, timesteps):

    N = len(matrix)
    S_A = np.zeros(timesteps)
    qubits = list(range(N))
    for i in range(timesteps):
        # print(matrix, 'initial:')
        q = random.sample(qubits, 1)[0]
        Tgate(matrix, q)
        del qubits[N-2:N]
        target = random.sample(qubits, 1)[0]
        control = random.sample([target, target + 1, target + 2], 1)[0]
        C3gate(matrix, target, control)
        submatrix = getCutStabilizers(matrix, round(N/2))
        S_A[i] = matrix_rank(submatrix) - N/2
    return S_A * math.log(2) / math.log(2**(N/2))
    # return matrix, S_A


def CNOT(q, matrix):
    N = len(matrix)
    xvec = matrix[:, :N]
    zvec = matrix[:, N: 2*N]

    xvec[q] += (xvec[q+1])
    xvec[q+1] = (xvec[q+1])
    zvec[q] = zvec[q]
    zvec[q+1] += zvec[q]
    matrix[:, :N] = xvec % 2
    matrix[:, N: 2*N] = zvec % 2

    return matrix


# def randomclifford(matrix, timesteps):
#     N = len(matrix)
#     S_A = -1 * round(N/2) * np.ones(timesteps)
#     qubits = list(range(N))
#     for i in range(timesteps):
#         j = random.sample(qubits, 1)[0]
#         choice = random.sample([1, 2], 1)[0]
#         if choice == 1:
#             Tgate(matrix, j)
#         elif choice == 2:
#             del qubits[N-1:N]
#             CNOT
#         S_A[i] += matrix_rank(matrix[:, N:])
#     return matrix, S_A


L = 100
tsteps = 20000

t = range(tsteps)


Entropy = {
    "Time": t,
    "Ent1": randomcircuit(checkmatrix(L), tsteps),
    "Ent2": randomcircuit(checkmatrix(L), tsteps),
    "Ent3": randomcircuit(checkmatrix(L), tsteps),
    "Ent4": randomcircuit(checkmatrix(L), tsteps),
    "Ent5": randomcircuit(checkmatrix(L), tsteps),
    "Ent6": randomcircuit(checkmatrix(L), tsteps),
    "Ent7": randomcircuit(checkmatrix(L), tsteps),
    "Ent8": randomcircuit(checkmatrix(L), tsteps),
    "Ent9": randomcircuit(checkmatrix(L), tsteps),
    "Ent10": randomcircuit(checkmatrix(L), tsteps),
}
df2 = pd.DataFrame(Entropy)
with pd.ExcelWriter('StabilizerCircuits100.xlsx') as writer:
    df2.to_excel(writer, sheet_name='Sheet_1')

print('done')


# output = randomcircuit(testmatrix, timesteps)[1]

# print(testmatrix)
# # print(matrix_rank(testmatrix[0:round(L/2), :]))
# output = randomcircuit(testmatrix, timesteps)[1]
# print("Page Value:", math.log(2**(L/2)))
# plt.plot(range(timesteps), output)
# plt.show()
