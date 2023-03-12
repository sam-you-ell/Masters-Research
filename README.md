# Quantum Scrambling of Local Operators

# Aims

# Methods

# Some Basic Theory

# Some Results


# Determination of the SO(m) matrices, 

Function to determine R comes from the following, 
1. Determine eigenvalues of $\alpha$
1. Block diagonalise $\alpha$, the matrix that forms our gate hamiltonian, 
2. Take the matrix, $W$ that block-diagonalises $\alpha$, and store. 
4. $R$ is then found by $R = W^{T} M W $. Where $M$ is a block diagonal matrix, formed from    the eigenvalues of $\alpha$.


# Evolving the System

To quantify the state of the system, such that the entanglement entropy can be calculated, we need
to determine the correlation matrix, $$M_{ij} = \langle \psi | U^T c_i c_j U |\psi \rangle$$, where $U$ is our
quantum circuit, and $\psi$ is taken to be the vacuum, denoted $|\underline{0}\rangle$. Then using
the result from DiVencenzo, we can recast this to, 
$$M_{ij} = \langle \psi | U^T c_i c_j U |\psi \rangle = \sum_{k l} R_{ik}R_{jl}\langle \underline{0} | c_i c_j  |
\underline{0} \rangle$$

# Correlation matrix Form

Correlation matrix, is given by

$$C_{m, n} = \frac{\sin(q_F (m-n))}{\pi (m-n)}$$

# Poster

Introduction

Stabilizer circuits
Free Fermion CIrcuits

Results
- Stabilizer entang entropy
- Free Fermion entang entropy


Conclusions


What we want to explain, 

What quantum scrambling
Combination of operator scrambling and entanglement growth






