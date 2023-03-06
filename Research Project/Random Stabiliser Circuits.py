import numpy as np
from scipy import sparse
import random
import sympy
from numpy.linalg import matrix_rank
import matplotlib.pyplot as plt


def checkmatrix(N):
    """
    Function sets up the check matrix, for the simplest case, will likely need to change so I can
    specify a string. 
    """
    return sparse.diags([0, 1, 0], [-1, 0, 1], shape=(N, 2*N)).toarray()
    # X = np.zeros((N, N))
    # Z = sparse.diags([0, 1, 0], [-1, 0, 1], shape=(N, N)).toarray()
    # return np.concatenate((X, Z), axis=1)


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
    zvec[q] = (zvec[q] + xvec[q+1] + zvec[q+1] + xvec[q+2] + zvec[q+2]) % 2
    xvec[q+1] = (xvec[q] + xvec[q+1]) % 2
    zvec[q+1] = (zvec[q+1] + xvec[q]) % 2
    xvec[q+2] = (xvec[q] + xvec[q+2]) % 2
    zvec[q+2] = (zvec[q+2] + xvec[q]) % 2

    matrix[:, :N] = xvec
    matrix[:, N: 2*N] = zvec

    if control == q:
        pass
    elif control == q+1:
        matrix[[q+1, q]] = matrix[[q, q+1]]
    elif control == q+2:
        matrix[[q+2, q]] = matrix[[q, q+2]]

    return matrix


def randomcircuit(matrix, timesteps):
    S_A = np.zeros(timesteps)
    N = len(matrix)
    p = 30
    for i in range(timesteps):
        # print(matrix, 'initial:')
        qubits = list(range(N))
        k = random.sample(qubits, 1)[0]
        Tgate(matrix, k)

        del qubits[N-2:N]

        j = random.sample(qubits, 1)[0]
        control = random.sample([j, j+1, j+2], 1)[0]
        C3gate(matrix, j, control)
        submatrix = matrix[:, 0:N]
        S_A[i] = matrix_rank(submatrix, tol=1)
    return matrix, S_A


testmatrix = checkmatrix(120)
print(testmatrix.shape)

timesteps = range(10000)
entr = randomcircuit(testmatrix, 10000)[1]

# print(randomcircuit(testmatrix, 2000)[1])

plt.plot(timesteps, entr)
plt.show()
