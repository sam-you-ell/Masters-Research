import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Clifford
import random
from qiskit.quantum_info.operators import PauliTable
from numpy.linalg import matrix_rank
import matplotlib.pyplot as plt


def checkmatrix(L, paulilist):
    matrix = np.zeros(shape=(L, 2*L))
    for i in range(len(paulilist)):
        for j in range(len(paulilist[i])-1):
            if paulilist[i][j+1] == 'I':
                matrix[i, j] = 0
                matrix[i, j+L] = 0
            elif paulilist[i][j+1] == 'X':
                matrix[i, j] = 1
                matrix[i, j+L] = 0
            elif paulilist[i][j+1] == 'Y':
                matrix[i, j] = 1
                matrix[i, j+L] = 1
            elif paulilist[i][j+1] == 'Z':
                matrix[i, j] = 0
                matrix[i, j+L] = 1
    return matrix


# N = 10
# timesteps = 1
# qc = QuantumCircuit(N)
# qubits = list(range(N))
# S_A = np.zeros(timesteps)
# for i in range(timesteps):
#     j = random.sample(qubits, 1)[0]
#     choice = random.sample([1, 2], 1)[0]
#     if choice == 1:
#         qc.h(j)
#     elif choice == 2:
#         del qubits[N-1:N]
#         q = random.sample(qubits, 1)[0]
#         qc.cx(0, q+1)
#     cliff = Clifford(qc)
#     mat = checkmatrix(N, cliff.to_labels(mode="S"))
#     S_A[i] = matrix_rank(mat[:, 0:N])
N = 120
timesteps = 10000


def bal_random(L, timesteps):
    S_A = np.zeros(timesteps)
    qubits = list(range(L))
    qc = QuantumCircuit(L)

    for i in range(timesteps):
        control, target, q = random.sample(qubits, 3)
        qc.h(q)
        qc.z(q)
        qc.cy(control, target)
        qc.cy(control, target)

        cliff = Clifford(qc)
        mat = checkmatrix(N, cliff.to_labels(mode="S"))
        S_A[i] = matrix_rank(mat[:, 0:N])
    return mat, S_A


output = bal_random(N, timesteps)[1]
# print(output)
plt.plot(range(timesteps), output)
plt.show()
