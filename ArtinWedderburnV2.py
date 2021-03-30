import numpy as np
import scipy.sparse as sparse


class SparseAlgebra:
    def random_vector(self):
        d = self.dimension
        return (np.random.rand(d) - 0.5 ) + 1j * (np.random.rand(d) - 0.5 )

    def commutator_matrix(self):
        d = self.dimension


        l = sparse.hstack(
            [lm.reshape(d * d, 1)
             for lm in self.left_multiplication_matrices])

        r = sparse.hstack(
            [rm.reshape(d * d, 1)
             for rm in self.right_multiplication_matrices])


        return l - r


    def multiplication_matrix(self,v, ms):
        # ms is a list of sparse matrices
        d = self.dimension
        accumulator = np.zeros((d,d), dtype=complex)

        for (c, m) in zip(v, ms):
            accumulator += c * m.todense()

        return accumulator

    def left_multiplication_matrix(self,v):
        return self.multiplication_matrix(v, self.left_multiplication_matrices)

    def right_multiplication_matrix(self,v):
        return self.multiplication_matrix(v, self.right_multiplication_matrices)



    def __init__(
            self,
            dimension,
            left_multiplication_matrices,
            right_multiplication_matrices,
            unit):

        # dimension of algebra
        # we denote the basis vectors as e_1, e_2, ...., e_dimension
        self.dimension = dimension

        # list of matrices
        # jth column of ith matrix is e_i e_j
        self.left_multiplication_matrices = left_multiplication_matrices

        # list of matrices
        # jth column of ith matrix is e_j e_i
        self.right_multiplication_matrices = right_multiplication_matrices

        # unit of algebra. Dense numpy array
        self.unit = unit


def load_sparse_algebra_from_file(path):
    with open(path,'r') as f:
        lines = f.readlines()

    unit = np.array([complex(x) for x in lines[0].split(' ')])
    dimension = len(unit)


    multiplication_lines = lines[1:]
    left_multiplication_data = [([],([],[])) for _ in range(dimension)]
    right_multiplication_data = [([],([],[])) for _ in range(dimension)]

    for line in multiplication_lines:
        parts = line.split(' ')
        i = int(parts[0])
        j = int(parts[1])
        k = int(parts[2])
        val_real = float(parts[3])
        val_imaginary = float(parts[4])
        val = complex(val_real, val_imaginary)


        left_multiplication_data[i][0].append(val)
        left_multiplication_data[i][1][0].append(k)
        left_multiplication_data[i][1][1].append(j)


        right_multiplication_data[j][0].append(val)
        right_multiplication_data[j][1][0].append(k)
        right_multiplication_data[j][1][1].append(i)

    left_multiplication_matrices = [
        sparse.coo_matrix(m, dtype=complex, shape=(dimension,dimension))
        for m in left_multiplication_data]

    right_multiplication_matrices = [
        sparse.coo_matrix(m, dtype=complex, shape=(dimension,dimension))
        for m in right_multiplication_data]

    return SparseAlgebra(
        dimension,
        left_multiplication_matrices,
        right_multiplication_matrices,
        unit)

