
import numpy as np
from scipy.sparse import csr_matrix, identity, kron
# Pauli-Spin Operators,


def pauli_spin(sigma: str):
    if sigma == 'X':
        return 1/2 * np.array(((0, 1), (1, 0)))

    elif sigma == 'Y':
        return 1/2 * np.array(((0, -1j), (1j, 0)))

    elif sigma == 'Z':
        return 1/2 * np.array(((1, 0), (0, -1)))

    elif sigma == 'I':
        return np.eye(2)


def spinops(L):
    X = 1/2 * csr_matrix([[0, 1], [1, 0]])
    Y = 1/2 * csr_matrix([[0, -1j], [1j, 0]])
    Z = 1/2 * csr_matrix([[1, 0], [0, -1]])
    I = identity(2)
    Sx = np.array(np.zeros((2*L, 2*L)))
    Sy = np.zeros((2*L, 2*L))
    Sz = np.zeros((2*L, 2*L))

    for j in range(L):
        Sx[j] = X
        # Sy[j] = Y
        # Sz[j] = Z
        for _ in range(j):
            Sx[j] = kron(I, Sx[j])
        for _ in range(j+1, L):
            Sx[j] = kron(Sx[j], I)

    return Sx


print(spinops(5))
