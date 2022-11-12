import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.providers.aer import QasmSimulator
from qiskit.visualization import plot_histogram

circ = QuantumCircuit(3)

circ.h(0)
circ.cx(0, 1)
circ.cx(0,2 )

print(circ.draw())

