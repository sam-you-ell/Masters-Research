import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.providers.aer import QasmSimulator
from qiskit.visualization import plot_histogram
from qiskit.circuit.library import XGate, YGate, ZGate

# circ = QuantumCircuit(3)

# circ.h(0)
# circ.cx(0, 1)
# circ.cx(0,2 )

# print(circ.draw())
#defining matrices
pX = XGate(); pY = YGate(); pZ = ZGate()
print(pX.to_matrix())
print(pY.to_matrix())
print(pZ.to_matrix())

# GHZcircuit = QuantumCircuit(3)

# GHZcircuit.h(0)
# GHZcircuit.cx(0, 1)
# GHZcircuit.cx(0, 2)
# print(GHZcircuit.draw())



