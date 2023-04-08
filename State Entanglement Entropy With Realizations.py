# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 12:02:16 2022

@author: lenovo
"""
import matplotlib.colors
import matplotlib as mpl
import numpy as np
import numpy.linalg as la
from qiskit import QuantumCircuit, assemble, Aer
from qiskit.visualization import plot_histogram
from qiskit import transpile
from qiskit.providers.aer import AerSimulator
import random
from qiskit.quantum_info import Clifford
import matplotlib.pyplot as plt


# =============================================================================
# =============================================================================
# # SIMULATING ENTANGLEMENT ENTROPY OF QUBIT SYSTEMS AVERAGED OVER REALIZATIONS
# =============================================================================
# =============================================================================

# Define Functions:
# =============================================================================
# Unitary time evolution, applying a random gate on random qubit(s)
# ============================================================================

def ZH(circuit, i, state=bool):
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


def SWAP(circuit, i, j, state=bool):
    """
    Two-qubit gate in state / operator (0 -> X, 1 -> Y) space:
    ----------
    Parameters
    ----------
    circuit : Qiskit Quantum Circuit.
    i,j : type=int, Qubits.
    state : type=bool, True for State Space, False for Operator Space

    """
    if state == True:
        # circuit.swap(i,j)
        # circuit.x(i)
        # circuit.y(j)
        # circuit.swap(i,j)
        circuit.y(i)
        circuit.x(j)
    if state == False:
        circuit.swap(i, j)
    return circuit


# =============================================================================
# C3 Gate
# =============================================================================
# State space, C3 = CX_21 CX_31 CZ_12 T^6_1 T^6_2
# Operator space, 0 -> X, 1 -> Y, C3 = CY_12 CY_13


def C3(circuit, i, j, k, state=bool):
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
        circuit.cx(j, i)
        circuit.cx(k, i)
        circuit.cz(i, j)
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
        circuit.cx(i, j)
        circuit.s(j)
        # CY_13
        circuit.sdg(k)
        circuit.cx(i, k)
        circuit.s(k)
    return circuit


def U(circuit, n, evolutions, text=False):
    """
    Minimum two-qubit circuit evolution:
    ----------
    Parameters
    ----------
    circuit : Qiskit Quantum Circuit.
    evolutions: type=int, Number of random gates in circuit.
    """
    if text == True:
        print("Evolving circuit with "+str(evolutions)+" random gates:")
    # Time evolutions
    steps = np.arange(0, evolutions, 1)
    for i in steps:
        # Indexing qubits
        qubit_indices = [m for m in range(n)]
        # Random 2 qubit selection
        i, j, k = random.sample(qubit_indices, 3)
        # Random gate
        rgate, = random.sample([1, 2, 3], 1)
        # print(rgate)
        if rgate == 1:
            ZH(circuit, i, state=False)
        elif rgate == 2:
            SWAP(circuit, i, j, state=False)
        elif rgate == 3:
            C3(circuit, i, j, k, state=False)
    return circuit

# =============================================================================
# Check Matrix function
# =============================================================================


def Check_Matrix(circuit, n, text=False):
    clifford_circuit = Clifford(circuit)
    Stabilizers = clifford_circuit.to_labels(mode="S")
    l = np.arange(0, n)
    index = np.arange(1, n+1)
    Paulis = ['I', 'X', 'Z', 'Y']
    checkmatrix = np.zeros((n, 2*n))
    if text == True:
        print()
        print("Check Matrix:")
    for j in l:
        checkx = np.zeros(len(index))
        checkz = np.zeros(len(index))
        for i in index:
            si = Stabilizers[j][i]
            if si == Paulis[1]:
                checkx[i-1] = 1
                checkmatrix[j][i-1] = 1
            if si == Paulis[2]:
                checkz[i-1] = 1
                checkmatrix[j][n+i-1] = 1
            if si == Paulis[3]:
                checkx[i-1] = 1
                checkz[i-1] = 1
                checkmatrix[j][i-1] = 1
                checkmatrix[j][n+i-1] = 1
        if text == True:
            print("g_"+str(j+1)+" = "+str(checkx)+'|'+str(checkz))
    return checkmatrix

# =============================================================================
# Submatrix of transposed check matrix for Region A, consisting of |A| qubits
# =============================================================================


def Submatrix_A(n, cmt, text=False):
    # Cut at qubit a, default is round(n/2)
    a = round(n/2)
    # modA = number of qubits in Region A
    modA = n - a
    # Submatrix is 2*|A| by n matrix
    submatrixA = np.zeros((2*modA, n))
    submatrixA[0:modA, :] = cmt[a:a+modA, :]
    submatrixA[modA:2*modA, :] = cmt[n+a:2*n, :]
    if text == True:
        print("Cut at qubit "+str(a) +
              "; Region A (right of cut) contains "+str(modA)+" qubits.")
    return submatrixA, modA
# =============================================================================
# Entanglement Entropy Calculation from submatrix and number of qubits
# =============================================================================


def S(submatrix, mod, text=False):
    I_A = la.matrix_rank(submatrix)
    s_A = I_A - mod
    if text == True:
        print("Entanglement entropy of Region A, S_A, = "+str(s_A))
    return s_A

# =============================================================================
# Entanglement Entropy Calculation function for set time
# =============================================================================


def Entanglement_Entropy(n, evolutions, text=bool):
    # Quantum circuit
    qc = QuantumCircuit(n, n)
    if text == True:
        print("Quantum circuit with "+str(n)+" qubits.")
        print("Evolving circuit with "+str(evolutions)+" random gates:")
    # Applying random unitary gates 'evolutions'# of times:
    U(qc, n, evolutions, text=False)
    # Getting check matrix for circuit:
    CM = Check_Matrix(qc, n, text=False)
    CMT = CM.T
    if text == True:
        print()
        print("Check Matrix (Transposed) for Stabilizers of System:")
        print(CMT)
        print()
    submatrixA, modA = Submatrix_A(n, CMT, text=text)
    if text == True:
        print()
        print("Submatrix of Region A:")
        print(submatrixA)
        print()
    S_A = S(submatrixA, modA, text=text)
    return S_A

# Demo Experiment:
# Number of qubits in quantum circuit (minimum of 2)
# n = 4
# Time evolutions (number of random gates applied)
# evolutions = 0
# Entanglement_Entropy(n,evolutions,text=True)

# =============================================================================
# Experiment over range of timesteps, calculating entanglement entropy at each
# =============================================================================


def Experiment_Over_Times(n, timesteps):
    # Array to store entanglements entropy values in
    S_As = np.zeros(len(timesteps))
    # Run simulation over range of times (specified in 'timesteps'):
    for i in range(len(timesteps)):
        S_As[i] = Entanglement_Entropy(n, timesteps[i], text=False)
    return S_As


# =============================================================================
# =============================================================================
# =============================================================================
# Full Simulation; averaged over realizations (default=1), plotting optional
# =============================================================================
# =============================================================================
# =============================================================================
# All simulated for 120 qubits
n = 120
# Specify number of time steps and total time
steps = 500
total_time = 1600

mpl.rcParams['text.usetex'] = False  # not really needed


def Full_Simulation(N, ts, T, Realizations=1, Plotting=True):
    """
    """
    # Make timesteps into list of integers up to total time
    t = np.linspace(0, T, ts)
    timesteps = []
    for i in range(len(t)):
        timesteps.append(int(t[i]))

    # Averaging over realizatons
    runs = Realizations
    # Array for many realizations
    S_A_Values = np.zeros((runs, len(timesteps)))
    for i in range(runs):
        S_As = Experiment_Over_Times(n, timesteps)
        S_A_Values[i, :] = S_As
    # Final array for realization-averaged entanglement entropies
    avS_As = np.zeros(len(timesteps))
    for i in range(len(timesteps)):
        avS_As[i] = np.mean(S_A_Values[:, i])

    if Plotting == True:
        # import matplotlib as mpl
        # mpl.rcParams['text.usetex'] = True
        # mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

        # plt.title(r"$f_{\mathrm{cor, r}}$")
        # plt.title(r"$f_{\mathrm{cor, r}}$")
        fig, ax = plt.subplots()
        ax.plot(timesteps, avS_As)

        # No. of time evolutions (H and CNOT)")
        # plt.title(r"$S_{A}(t)$ $vs$ no. of time evolutions (applying H and CNOT)", fontstyle='italic')
        plt.xlabel(r"  Time", rotation=0, loc='center')
        plt.ylabel(r" $S_{A}(t)$         ", rotation=0, loc='center')
        # ax.grid()
        ax.grid(which='major', color='dimgray', linewidth=0.1)
        ax.grid(which='minor', color='lightgray', linewidth=0.1)
        ax.minorticks_on()
        plt.xlim([0, 1600])
        plt.ylim([0, 62])
        # removing top and right borders
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.tight_layout()
        plt.show()

    return avS_As, timesteps


Average_S_As, Times = Full_Simulation(n, steps, total_time, Realizations=1)
