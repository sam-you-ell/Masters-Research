import numpy as np
from scipy.sparse import csr_matrix, identity, kron, lil_matrix
from scipy.stats import unitary_group
import scipy
import math

U = unitary_group.rvs(2)
np.set_printoptions(precision=3, suppress=True)


def basis(L: int, flip: bool):
    """
    Creates all possible configurations for a given spin-chain size L. Outputs an (2^L x L) matrix
    of configurations.
    """
    space = np.zeros((2**L, L)
                     )  # initialise vector of basis vectors, 2^L x L space.
    # just supresses scientific notation when printing, no effect on floating point errors
    np.set_printoptions(suppress=True)
    binary_numbs = []
    for i in range(2**L):  # main method. converts the index number, i into binary form using an f string. then stores the binary string in our space
        binary_numbs.append(f'{i:0{L}b}')
        for j in range(L):
            space[i, j] = binary_numbs[i][j]
    # can flip if needed. seems to depend on how you order your basis states , this needs clarification***
    return np.fliplr(space) if flip else space


# def occupation_num()


def creation(L: int, site: int):
    """
    Function that creates the creation operator from second quantisation (in matrix form). The
    functions takes a chain length:L and the index of the operator:site and returns the operator
    matrix.

    Method explained in code.
    """
    matrix = lil_matrix((2**L, 2**L))
    input_basis = basis(L, False)
    output_basis = np.zeros((2**L, L))
    c = np.zeros(L)
    c[site] = 1

    for i in range(2**L):  # loops over entire space
        # add one to the element at specified site. If greater than 1, unphysical
        output_basis[i] = input_basis[i] + c
        # if any(output_basis[i, j] > 1 for j in range(L)):
        #     output_basis[i] = np.empty(L).fill(np.nan) # filters out unphysical states, not sure
        #     if needed.
        # if state is physical, we check if it was in the original basis, and find the index in both bases.
        for j in range(len(input_basis)):
            if (output_basis[i] == input_basis[j]).all():
                # print(i, j)
                # fills element in matrix with the phase factor, as specified in the formalism.
                matrix[j, i] = (-1)**(np.sum(input_basis[i]))
    return matrix.tocsr()


# print(creation(3, 0))

# checking anti-commutation relations


def anti_commutator(matrix1, matrix2, majorana):
    """
    function to check if creation/annihilation operators obey the correct (fermionic) canonical
    anti-commutation rules, i.e {a_i, a_j} = {a_i^dag, a_j^dag} = 0 and {a_i, a_j^dag} = delta_ij
    Function takes two matrices and returns the matrix after applying the anti-commutator operation.
    """
    if majorana is True:
        return (np.matmul(matrix1.toarray(), matrix2.toarray()) + np.matmul(matrix2.toarray(), matrix1.toarray()))
    else:
        return (np.matmul(matrix1.toarray(), matrix2.toarray()) - np.matmul(matrix2.toarray(), matrix1.toarray())) % 2
    # return ((matrix1.multiply(matrix2)) + x * (matrix2.multiply(matrix1))).toarray()


def majorana(L, j, even):
    """
    Defining the majorana operators for even and odd sectors. Even sector is indexed by 2j while odd
    sector is indexed by 2j+1
    """
    m = 2*L
    c = lil_matrix((2**m, 2**m))
    return creation(L, j) + creation(L, j).transpose() if even is True else -1j*(creation(L, j).transpose() - creation(L, j))


def gen_gate_hamil(L):
    b = np.zeros((L, L))
    b[1, 1] = -1
    b[0, 0] = 1
    H_g = np.zeros((2**L, 2**L))
    for i in range(L):
        for j in range(L):
            H_g += b[i, j] * creation(L, i) @ creation(L, j).transpose()
    return H_g

# c1 c2 c3 c4
# c0 c1 c2 c3


# def UH1():
#     H = -(math.pi/4) * (-majorana(2, 0, True) @ majorana(2, 1, False) + majorana(2, 0, False) @ majorana(2, 1,
#                                                                                                          True) + majorana(2, 0, True) @ majorana(2, 0, False) + majorana(2, 1, True) @ majorana(2, 1, False))
#     return H


def ent():
    H = -(math.pi/8) * (majorana(2, 0, True)@majorana(2, 0, False))
    return scipy.linalg.expm(1j*H)


def correlation(L):
    m = 2*L
    M = lil_matrix((L, L))
    vac = np.zeros(shape=(2**L, 1))
    vac[1] = 1

    for i in range(L):
        for j in range(L):

            M[i, j] = 1j*(vac.transpose() @ majorana(L, i, True) @ majorana(L, j, False) @ vac + vac.transpose() @ majorana(L, i, True) @ majorana(L, j, False) @ vac +
                          vac.transpose() @ majorana(L, i, False) @ majorana(L, j, True) @ vac +
                          vac.transpose() @ majorana(L, i, False) @ majorana(L, j, False) @ vac - vac.transpose() @ majorana(L, j, True) @ majorana(L, i, False) @ vac + vac.transpose() @ majorana(L, j, True) @ majorana(L, i, False) @ vac +
                          vac.transpose() @ majorana(L, j, False) @ majorana(L, i, True) @ vac +
                          vac.transpose() @ majorana(L, j, False) @ majorana(L, i, False) @ vac
                          )
    return M.todense()


N = 4


print(correlation(N))


# creation_0 = creation(N, 0)
# creation_1 = creation(N, 1)
# c_0 = majorana(2, 0, True).toarray()
# print(c_0)
# print(UH1().toarray())
# print(ent().conj().transpose() @ c_0 @ ent())
# # print(gen_gate_hamil(2))
# # L = 12
# # print(majorana(L, 1, True))
# # print(-1j*np.matmul(majorana(2, 1, False).toarray(), majorana(2, 1, True).toarray()))
# # define some hamiltonians,
# ################
# # # Checking Commutation Relations
# N = 3
# # a_dag_1 = creation(N, 0)
# # a_dag_2 = creation(N, 1)
# # a_1 = creation(N, 0).transpose()
# # # print(a_dag_1.toarray())
# # c_2 = majorana(N, 1, True)

# # c_3 = majorana(N, 1, False)
# # print(a_dag_1.toarray())

# vac = np.zeros(shape=(2**N, 1))
# vac[1] = 1


# print(vac)
# a_dag_2 = creation(N, 1).toarray()
# fock = np.matmul(a_dag_2, vac)
# print(fock)
# print(a_dag_2)
# print(anti_commutator(c_2, c_3))

# VEV = vac.transpose() @ (anti_commutator(creation(N, 0).transpose(),
#                                          creation(N, 0).transpose(), False)) @ vac
# print(VEV)
# # print(a_1.toarray())
# # print(a_dag_1.toarray())
# # doesnt give right anti commutation relation. I'm not sure why
# print(anti_commutator(a_dag_2, a_1, False))


#################

# x = creation(N, 0).toarray() @ creation(N, 0).toarray().transpose()
# y = 1/2 * (np.eye(2**N) + (1j * majorana(N, 0, True).toarray()
#                            @ majorana(N, 0, False).toarray()))

# print(x)
# print(y)

# correlation_mat = vac.transpose() @ majorana(N, )
# print(Data[str(basis(m)[2])])
# idx = np.arange(2**L)
# Data = {str(input_basis[i]): idx[i] for i in range(2**L)}
# num_bits = 4  # L
# num = 2  # Input number - i in my case
# print(f'{num:0{num_bits}b}')  # output


# print(unitary_group.rvs(2))
