import numpy as np
from scipy import linalg

# class Pauli:

#     def __init__(self, X, Y, Z):
#         self.X = X
#         self.Y = Y
#         self.Z = Z

#     X = np.array(([0j, 1 + 0j], [1+0j, 0j]))
#     Y = np.array(([0j, (-1)*1j], [1j, 0j]))
#     Z = np.array(([1+0j, 0j], [0j, -1+0j]))


# class Unitary(Pauli):
#     def __init__(self, H, Toff):
#         self.H = H #Hadamard
#         self.Toff = Toff#Toffoli Phase

#     H = (1 / np.sqrt(2)) * np.array(([1, 1],[1, -1]))
#     Toff = np.array(([1, 0], [0, (1 + 1j)/np.sqrt(2)]))

#     def transform(self):
#         return np.matmul(self.conj().T, np.matmul(Pauli.X, self))

# res = Unitary.transform(Unitary.H)
# print(res)

# def unitarytransform(pauli, unitary):
#     pauli = np.zeros(2, 2)
#     if unitary == H
#     transform =  np.matmul(unitary, np.matmul(pauli, np.matrix.transpose(unitary) ))
# print(trans)

##############################
from scipy.sparse import diags


def checkmatrix(N):
    return diags([0, 1, 0], [-1, 0, 1], shape=(N, 2*N)).toarray()


print(checkmatrix(3))


def had(chm, N, i, j):
    # chm = np.zeros((N, 2*N))
    hada = np.zeros(2*N)
    hada[j] = 1
    hada[j + N] = 1
    chm[i] = (chm[i] + hada) % 2
    return chm


print(had(checkmatrix(4), 4, 1, 1))
