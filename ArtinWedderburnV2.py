import numpy as np
import scipy.sparse as sparse
from scipy.linalg import eigh

class ArtinWedderburn:
    def compute_center(self):
        algebra = self.sparse_algebra
        d = algebra.dimension

        cm = algebra.commutator_matrix()
        cm_transpose = cm.transpose()
        cm_conj = cm_transpose.conj()
        cm_conj_times_cm = cm_conj * cm
        s, v = eigh(cm_conj_times_cm.todense())
        center_dimension = len(np.where( np.abs(s) < self.threshold)[0])
        center_inclusion = v[:,:center_dimension]
        center_projection = np.transpose(np.conj(center_inclusion))

        center_multiplication = [[ np.dot(
            center_projection,
            algebra.multiply(v[:,i], v[:,j]))
           for i in range(center_dimension)]
         for j in range(center_dimension)]


    def log(self,*args, **kwargs):
        if self.logging:
            print(*args, **kwargs)

    def __init__(self, sparse_algebra, threshold = 1.0e-5, logging = False):
        self.sparse_algebra = sparse_algebra
        self.threshold = threshold
        self.logging = logging
        self.total_defect = sparse_algebra.algebra_defect()

        self.log("algebra defect:", self.total_defect)
        self.log("")

        self.log("computing center...")
        self.compute_center()
        self.log("")



class SparseAlgebra:
    def associative_defect(self):
        x = self.random_vector()
        y = self.random_vector()
        z = self.random_vector()
        l = self.multiply(self.multiply(x,y), z)
        r = self.multiply(x, self.multiply(y,z))
        return np.sum(np.abs(l - r))

    def left_identity_defect(self):
        lm = self.left_multiplication_matrix(self.unit)
        return np.sum(np.abs(lm - np.identity(self.dimension)))

    def right_identity_defect(self):
        rm = self.right_multiplication_matrix(self.unit)
        return np.sum(np.abs(rm - np.identity(self.dimension)))

    def algebra_defect(self):
        return sum([self.associative_defect(),
                    self.left_identity_defect(),
                    self.right_identity_defect()])

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

    def helper(self,v, ms):
        # ms is a list of sparse matrices
        d = self.dimension
        accumulator = np.zeros((d,d), dtype=complex)

        for (c, m) in zip(v, ms):
            accumulator += c * m.todense()

        return accumulator

    def left_multiplication_matrix(self,v):
        return self.helper(v, self.left_multiplication_matrices)

    def right_multiplication_matrix(self,v):
        return self.helper(v, self.right_multiplication_matrices)

    def multiply(self, x, y):
        return self.multiplication_matrix * (np.kron(x,y)).flatten()

    def __init__(
            self,
            dimension,
            left_multiplication_matrices,
            right_multiplication_matrices,
            multiplication_matrix,
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


        self.multiplication_matrix = multiplication_matrix

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
    multiplication_data = ([],([],[]))

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

        multiplication_data[0].append(val)
        multiplication_data[1][0].append(k)
        multiplication_data[1][1].append(dimension * i + j)

    left_multiplication_matrices = [
        sparse.coo_matrix(m, dtype=complex, shape=(dimension,dimension))
        for m in left_multiplication_data]

    right_multiplication_matrices = [
        sparse.coo_matrix(m, dtype=complex, shape=(dimension,dimension))
        for m in right_multiplication_data]

    multiplication_matrix = sparse.coo_matrix(
        multiplication_data,
        dtype=complex,
        shape=(dimension, dimension * dimension))

    return SparseAlgebra(
        dimension,
        left_multiplication_matrices,
        right_multiplication_matrices,
        multiplication_matrix,
        unit)

