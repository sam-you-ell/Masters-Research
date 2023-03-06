import numpy as np
import random
import math

# ---------------------------------------------------------------------------------------------
# Can minus signs on array capture all info necessary for fock state to work only with arrays?
# ---------------------------------------------------------------------------------------------


def Vac(n):
    vac = np.zeros(n)
    return vac

# Fock states are arrays, with all values carrying the overall sign.


def Fock_State(occupations):
    """Returns the (unsigned) Fock state as a decorative string when given the occupations indices.
    occupations : type numpy.ndarray"""
    if type(occupations) != np.ndarray:
        raise Exception("occupations must be of type numpy.ndarray")
    ket = "|"
    for i in range(len(occupations)):
        ket = ket + str(int(abs(occupations[i])))
        if len(ket) < 2*len(occupations):
            ket = ket + ","
    sign = math.copysign(1, np.amin(occupations))
    s = '+'
    if sign == -1:
        s = '-'
    ket = s + ket + ">"
    return ket


# Define Creation and Annihiliation operators, which incorporate anti-commutation relations with the sign flips

def Creation(site, occupations):
    """Fermionic creation operator, acts as 'a' on specified site in given Fock state array.
    Returns final Fock state as array.
    site : type int
    occupations : type numpy.ndarray"""
    if type(site) != int:
        raise Exception("site must be of type int")
    elif site >= len(occupations):
        raise Exception("site must be < n")
    elif type(occupations) != np.ndarray and type(occupations) != float:
        raise Exception("occupations must be of type numpy.ndarray")
    out = 0.0
    if type(occupations) == np.ndarray:
        out = np.empty(len(occupations))
        for i in range(len(occupations)):
            out[i] = occupations[i]
        val = abs(out[site])
        sum = 0
        for i in range(site):
            sum += int(out[i])
        if val == 0:
            out[site] = math.copysign(1.0, np.amin(occupations))
            out = (-1)**sum * out
        elif val == 1:
            out = 0.0
    return out


def Annihiliation(site, occupations):
    """Fermionic annihiliation operator, acts as 'a^(dagger)' on specified site in given Fock array.
    Returns final Fock state as array.
    site : type int
    occupations : type numpy.ndarray"""
    if type(site) != int:
        raise Exception("site must be of type int")
    elif site >= len(occupations):
        raise Exception("site must be < n")
    elif type(occupations) != np.ndarray and type(occupations) != float:
        raise Exception("occupations must be of type numpy.ndarray")
    out = 0.0
    if type(occupations) == np.ndarray:
        out = np.empty(len(occupations))
        for i in range(len(occupations)):
            out[i] = occupations[i]
        val = abs(out[site])
        sum = 0
        for i in range(site):
            sum += int(out[i])
        if val == 1:
            out[site] = math.copysign(0.0, np.amin(occupations))
            out = (-1)**sum * out
        elif val == 0:
            out = 0.0
    return out

# Function to create fermions on m random sites in a Fock state, usually the Vaccuum state.


def Random_Creation(m, occupations):
    """Creates fermions at m random sites in given Fock array, returns final Fock array.
    m : type int, m <= n
    occupations : type numpy.ndarray"""
    if type(occupations) != np.ndarray:
        raise Exception("fockstate must be of type numpy.ndarray")
    elif type(m) != int:
        raise Exception("m must be of type int")
    elif m > len(occupations):
        raise Exception("site must be <= n")
    indices = [i for i in range(len(occupations))]
    sites = random.sample(indices, m)
    for i in range(m):
        occupations = Creation(sites[i], occupations)
    return occupations

# # Function to create fermions on sites specified by list argument


def Specify_State(n, occupied):
    """Returns a Fock State from the Vacuum with fermions at specified sites in given Fock state.
    Specified indices of occupied sites with argument 'occupied', list containing positions.
    All positions specified in list must be < n.
    occupied : type list"""
    for i in occupied:
        if i >= n:
            raise Exception("specified sites must be <= n")
    F = np.zeros(n)
    for i in range(n):
        if i in occupied:
            F = Creation(i, F)
    return F


# Inner product function

def Inner(fockstate1, fockstate2):
    delta = None
    if abs(fockstate1) == abs(fockstate2):
        delta = 1.0 * math.copysign(1, np.amin(fockstate1)) * \
            math.copysign(1, np.amin(fockstate2))
    else:
        delta = 0.0
    return delta


# Binomial Coefficients through Pascal's Triangle
def Pascal_Triangle(N, verbose=False):
    coeff = []
    for i in range(N+1):
        vals = []
        for j in range(0, i+1):
            vals.append(math.comb(i, j))
        if verbose:
            print(vals)
        coeff.append(vals)
    return coeff


# Function to make permutations of particles in state - same as flipped binary representations of 0 to m
def Permutations(n, m):
    occ = np.zeros((m, n))
    bit = []
    for i in range(m):
        bit.append(bin(i)[2:])
        while len(bit[i]) < n:
            bit[i] = '0'+bit[i]
        for j in range(n):
            occ[i, j] = bit[i][j]
    occ = np.fliplr(occ)
    return occ


n = 4
Particle_no = np.arange(0, n+1)
Vacuum = np.zeros(n)
Pascal_row = Pascal_Triangle(n)[n]
m = sum(Pascal_row)
Possible_Occupations = Permutations(n, m)
# # print(m)
# print(Possible_Occupations)
Coefficients = np.arange(m)
Data = {}
for i in range(m):
    Data[str(Possible_Occupations[i])] = Coefficients[i]
print(Data)


# b = [2,0]
# Fock = abs(Specify_State(n, b))
# print(Fock)
# print("^ Represents: "+Fock_State(Fock))
# print(Data[str(abs(Fock))])
# # print()

# H = np.zeros((n,n))
# for i in range(n-1):
#     H[i,i] = 1
#     # H[i,i+1] = H[i+1,i] = 1
# # print(H)

# # print(H @ Fock)

# # print()

# # Redefining Creation and Annihiliation for easier reading
# State = None
# def a(i,F=State):
#     F = Annihiliation(i,F)
#     return F
# def a_d(i,F=State):
#     F = Creation(i,F)
#     return F


# # See possible outcomes of a hopping Hamiltonian
# # for j in Possible_Occupations:
# #     Fock[:] = Possible_Occupations[:,j]
# # i = 0
# # while i < n-1:
# #     if Fock[i] == 0.0 and Fock[i+1] == 1.0:
# #         F = abs(a(i+1,a_d(i,Fock)))
# #         print(F)
# #         # print(Data[str(abs(F))])
# #     if Fock[i+1] == 0.0 and Fock[i] == 1.0:
# #         F = abs(a(i,a_d(i+1,Fock)))
# #         print(F)
# #         # print(Data[str(abs(F))])
# #     i += 1
# print(Possible_Occupations[0])

# for j in range(len(Possible_Occupations)):
#     # print(Possible_Occupations[j])
#     Fock = Possible_Occupations[j]
#     # print(Fock)
#     print()
#     i = 0
#     # while i < n:
#     if Fock[i] == 0.0:
#         F = abs(a_d(i,Fock))
#         print(F)
#         print(Data[str(abs(F))])
#         # i += 1


# print()


# # def a_d_list(Fock):
# #     a_ds = []
# #     for i in range(len(Fock)):
# #         if abs(Fock[i]) == 1:
# #             a_ds.append(i)
# #     return a_ds
# # print(a_d_list(Fock))
