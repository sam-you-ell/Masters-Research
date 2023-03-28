# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 15:35:48 2022

@author: lenovo
"""

# =============================================================================
# BL21 Super-Clifford Stabilizer
# =============================================================================

import numpy as np
from qiskit import QuantumCircuit, assemble, Aer
from qiskit.visualization import plot_histogram

# =============================================================================
# Quantum Circuits
# =============================================================================

# Defining functions that replicates the gates from BL21 using Qiskit
# i.e. Super-Clifford Operators


# =============================================================================
# Z.H Gate
# =============================================================================
# Z.H = Tdg X T = (1 / np.sqrt(2)) * (X - Y)
# Z.H = Tdg Y T = (1 / np.sqrt(2)) * (X + Y)
# 0 -> X, 1 -> Y.


def ZH(circuit,i,state=bool):
    """
    Single-qubit gate in state / operator (0 -> X, 1 -> Y) space:
    ----------
    Parameters
    ----------
    circuit : Qiskit Quantum Circuit.
    i : type=int, Qubit.
    state : type=bool, True for State Space, False for Operator Space

    """
    if state == True:
        # Action of T in state space
        circuit.t(i)
    elif state == False:
        # Heisenberg conjugation of operators by T gives Z.H
        circuit.h(i)
        circuit.z(i)
    return circuit



# =============================================================================
# SWAP
# =============================================================================
# SWAP = SWAPdg X_1 Y_2 SWAP = Y_1 X_2
# 0 -> X, 1 -> Y.


def SWAP(circuit,i,j,state=bool):
    """
    Two-qubit gate in state / operator (0 -> X, 1 -> Y) space:
    ----------
    Parameters
    ----------
    circuit : Qiskit Quantum Circuit.
    i,j : type=int, Qubits.
    state : type=bool, True for State Space, False for Operator Space

    """
    if state==True:
        # circuit.swap(i,j)
        # circuit.x(i)
        # circuit.y(j)
        # circuit.swap(i,j)
        circuit.y(i)
        circuit.x(j)        
    if state==False:
        circuit.swap(i,j)      
    return circuit


# =============================================================================
# C3 Gate
# =============================================================================
# State space, C3 = CX_21 CX_31 CZ_12 T^6_1 T^6_2
# Operator space, 0 -> X, 1 -> Y, C3 = CY_12 CY_13


def C3(circuit,i,j,k,state=bool):
    """
    Three-qubit gate in state / operator (0 -> X, 1 -> Y) space:
    ----------
    Parameters
    ----------
    circuit : Qiskit Quantum Circuit.
    i,j,k : type=int, Qubits.
    state : type=bool, True for State Space, False for Operator Space
    """
    # Action in state space
    if state == True:
        circuit.cx(j,i)
        circuit.cx(k,i)
        circuit.cz(i,j)
        m = 0
        while m < 6:
            circuit.t(i)
            circuit.t(j)
            m += 1
    # Action in operator space, taking 0 -> X, 1 -> Y
    elif state == False:
        # A controlled Y
        # circuit.sdg(t)
        # circuit.cx(c,t)
        # circuit.s(t)
        
        # CY_12
        circuit.sdg(j)
        circuit.cx(i,j)
        circuit.s(j)        
        # CY_13
        circuit.sdg(k)
        circuit.cx(i,k)
        circuit.s(k)             
    return circuit


# =============================================================================
# Testing Gates
# =============================================================================

from qiskit import transpile
from qiskit.providers.aer import AerSimulator
import random


# # Quantum circuit with n qubits
# n = 9
# qc = QuantumCircuit(n,n)

# # Apply some random CNOT and S gates over m qubits
# qubit_indices = [i for i in range(n)]
# m = 1
# for i in range(m):
#     i, j, k, zh, swapi, swapj = random.sample(qubit_indices,6)
#     print(i,j,k)
#     print(zh)
#     ZH(qc,zh,state=False)
#     SWAP(qc,swapi,swapj,state=False)
#     C3(qc,i,j,k,state=False)


# qc.measure(range(n),range(n))

# ### Circuit Diagram
# qc.draw()
# print(qc)
# ### Simulations over circuit
# sim = Aer.get_backend('aer_simulator')
# qobj = assemble(qc)
# result = sim.run(qobj).result()
# counts = result.get_counts()
# print(counts)


# # Starting from vacuum string, the only super-Clifford operator that generates
# # any non-zero result is Z.H, as it takes to a super-position of X and Y



# =============================================================================
# Unitary time-evolutions from initial zero-string
# =============================================================================

# Quantum circuit with n qubits (multiple of 3, minimum is 3)
n = 40
qc = QuantumCircuit(n,n)


# Unitary time evolution, applying a random gate on random qubit(s)
def U(circuit,evolutions):
    """
    Minimum three-qubit circuit evolution in operator (0 -> X, 1 -> Y) space:
    ----------
    Parameters
    ----------
    circuit : Qiskit Quantum Circuit.
    evolutions: type=int, Number of random gates in circuit.
    """
    print("Evolving circuit with "+str(evolutions)+" random gates:")
    # Time evolutions
    steps=np.arange(0,evolutions,1)
    for i in steps:
        # Indexing qubits
        qubit_indices = [m for m in range(n)]
        # Random 3 qubit selection
        i, j, k = random.sample(qubit_indices,3)
        # print(i,j,k)
        # Random gate
        rgate, = random.sample([1,2,3],1)
        # print(rgate)
        if rgate == 1:
            ZH(circuit,i,state=False)
            print('Z.H on qubit '+str(i))           
        elif rgate == 2:
            SWAP(circuit,i,j,state=False)
            print('SWAP between qubits '+str(i)+' and '+str(j))           
        elif rgate == 3:
            C3(circuit,i,j,k,state=False)
            print('C3 on qubits '+str(i)+', '+str(j)+' and '+str(k))
    return circuit         


U(qc,8)

from qiskit.quantum_info import Clifford
cliff = Clifford(qc)

qc.measure(range(n),range(n))



### Circuit Diagram
qc.draw()
# print(qc)
### Simulations over circuit
sim = Aer.get_backend('aer_simulator')
qobj = assemble(qc)
result = sim.run(qobj).result()
counts = result.get_counts()
print(counts)


# Starting from vacuum string, the only super-Clifford operator that generates
# any non-zero result is Z.H, as it takes to a super-position of X and Y


print(cliff.to_labels(mode="S"))







