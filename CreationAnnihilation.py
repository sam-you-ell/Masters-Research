import numpy as np
from scipy.sparse import csr_matrix, identity, kron, lil_matrix


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


def anti_commutator(matrix1, matrix2):
    """
    function to check if creation/annihilation operators obey the correct (fermionic) canonical
    anti-commutation rules, i.e {a_i, a_j} = {a_i^dag, a_j^dag} = 0 and {a_i, a_j^dag} = delta_ij
    Function takes two matrices and returns the matrix after applying the anti-commutator operation.
    """

    return (np.matmul(matrix1.toarray(), matrix2.toarray()) + np.matmul(matrix2.toarray(), matrix1.toarray())) % 2
    # return ((matrix1.multiply(matrix2)) + x * (matrix2.multiply(matrix1))).toarray()


N = 2
a_dag_1 = creation(N, 0)
a_dag_2 = creation(N, 1)
a_1 = creation(N, 0).transpose()
print(a_dag_1.toarray())
# print(a_1.toarray())
# print(a_dag_1.toarray())
# print(commutator(a_1, a_dag_1))


#################

# print(Data[str(basis(m)[2])])
# idx = np.arange(2**L)
# Data = {str(input_basis[i]): idx[i] for i in range(2**L)}
# num_bits = 4  # L
# num = 2  # Input number - i in my case
# print(f'{num:0{num_bits}b}')  # output
