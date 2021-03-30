import numpy as np
import scipy.sparse as sparse
from scipy.linalg import eigh, eig

# takes an array of complex numbers and removes duplicates upto threshold
def fuzzy_filter(array, threshold):
    result = []
    for x in array:
        already_seen = False
        for v in result:
            if abs(x-v) < threshold:
                already_seen = True
                break

        if not already_seen:
            result.append(x)

    return result


def format_error(x):
    return str.format('{0:1.0e}' ,x)

class ArtinWedderburn:
    def compute_center(self):
        algebra = self.algebra
        d = algebra.dimension

        cm = algebra.commutator_matrix()
        cm_transpose = cm.transpose()
        cm_conj = cm_transpose.conj()
        cm_conj_times_cm = cm_conj * cm
        s, v = eigh(cm_conj_times_cm.todense())
        center_dimension = len(np.where( np.abs(s) < self.threshold)[0])
        center_inclusion = v[:,:center_dimension]
        center_projection = np.transpose(np.conj(center_inclusion))

        center_multiplication = np.array([[ np.dot(
            center_projection,
            algebra.multiply(v[:,i], v[:,j]))
           for i in range(center_dimension)]
         for j in range(center_dimension)])

        center = Algebra(
            center_dimension,
            center_multiplication,
            np.dot(center_projection, algebra.unit))

        self.center = center
        self.center_inclusion = center_inclusion
        center_defect = center.algebra_defect() + center.commutative_defect()
        self.log("center defect:", format_error(center_defect))
        self.total_defect += center_defect


    def compute_unscaled_central_idempotents(self):
        algebra = self.algebra
        center = self.center
        center_inclusion = self.center_inclusion
        v = center.random_vector()
        lv = center.left_multiplication_matrix(v)
        eigenvectors = eig(lv)[1]
        unscaled_central_idempotents = np.dot(center_inclusion, eigenvectors)

        unscaled_central_idempotent_defect = 0.0
        for i in range(center.dimension):
            for j in range(center.dimension):
                if i != j:
                    unscaled_central_idempotent_defect += np.abs(np.sum(algebra.multiply(
                        unscaled_central_idempotents[:,i],
                        unscaled_central_idempotents[:,j])))

        self.log("unscaled central idempotent defect:",
                 format_error(unscaled_central_idempotent_defect))
        self.total_defect += unscaled_central_idempotent_defect

        self.unscaled_central_idempotents = unscaled_central_idempotents


    def compute_block(self, idempotent_index):
        algebra = self.algebra
        idempotent = self.unscaled_central_idempotents[:,idempotent_index]
        left_multiplication = algebra.left_multiplication_matrix(idempotent)


        # u is the change of basis from new basis to old basis
        u,s,v =  svd(left_multiplication)

        # u is unitary, so inverse is conjugate transpose
        # u_inverse is change of basis from old basis to new basis
        u_inverse = np.transpose(np.conj(u))

        # the basis for the block is the columns of u whose singular value
        # is above the threshold

        counter = 0
        for i in range(len(s)):
            if abs(s[i]) > self.threshold:
                counter += 1

        block_inclusion = u[:,:counter]
        self.block_inclusions[idempotent_index] = block_inclusion


        block_dimension = block_inclusion.shape[1]
        m = algebra.multiplication
        d = algebra.dimension

        change_codomain_basis = np.tensordot(m,u_inverse, (2,1))
        change_left_basis = np.tensordot(u,change_codomain_basis,(0,0))
        m_rotated_anti = np.tensordot(u,change_left_basis,(0,1))
        m_rotated = np.transpose(m_rotated_anti, (1,0,2))

        m_defect =  np.sum(np.abs(m_rotated[0:block_dimension,
                                            0:block_dimension,
                                            block_dimension:d]))

        self.log("SVD multiplication defect:", format_error(m_defect))
        self.total_defect += m_defect

        unit_rotated = np.tensordot(idempotent, u_inverse, (0,1))

        unit_defect = np.sum(np.abs(unit_rotated[block_dimension:d]))
        self.log("SVD unit defect:", format_error(unit_defect))
        self.total_defect += unit_defect



    def log(self,*args, **kwargs):
        if self.logging:
            print(*args, **kwargs)

    def __init__(self, sparse_algebra, threshold = 1.0e-5, logging = False):
        self.algebra = sparse_algebra
        self.threshold = threshold
        self.logging = logging
        self.total_defect = sparse_algebra.algebra_defect()

        self.log("algebra defect:", self.total_defect)
        self.log("")

        self.log("computing center...")
        self.compute_center()
        self.log("")

        self.log("computing unscaled central idempotents...")
        self.compute_unscaled_central_idempotents()
        self.log("")

        self.blocks = {}
        self.central_idempotents = {}

        # the block inclusions are unitary, so to project onto a block
        # take the conjugate transpose of the inclusion
        self.block_inclusions = {}
        self.log("computing blocks...")
        self.log("")
        for i in range(self.center.dimension):
            self.log("computing block ", i)
            self.compute_block(i)
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
                    self.associative_defect(),
                    self.associative_defect(),
                    self.left_identity_defect(),
                    self.right_identity_defect()])

    def irrep_defect_multiplication(self,irrep):
        x = self.random_vector()
        y = self.random_vector()
        xy = self.multiplication(x,y)
        l = np.tensordot(xy, irrep, (0,0))
        r_uneval = np.transpose(np.tensordot(irrep, irrep, (2,1)),(2,0,1,3))
        r = np.tensordor(y, np.tensordot(x,r_uneval,(0,0)))
        result = np.sum(np.abs(l - r))
        return result


    def irrep_defect_identity(self,irrep):
        return np.sum(np.abs(np.tensordot(self.unit, irrep, (0,0)) - np.identity(irrep.shape[1])))

    def irrep_defect(self,irrep):
        return sum([ self.irrep_defect_multiplication(irrep),
                     self.irrep_defect_multiplication(irrep),
                     self.irrep_defect_multiplication(irrep),
                     self.irrep_defect_identity(irrep)])

    def multiply(self, x, y):
        return self.multiplication_matrix * (np.kron(x,y)).flatten()

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

    def mult_helper(self,v, ms):
        # ms is a list of sparse matrices
        d = self.dimension
        accumulator = np.zeros((d,d), dtype=complex)

        for (c, m) in zip(v, ms):
            accumulator += c * m.todense()

        return accumulator

    def left_multiplication_matrix(self,v):
        return self.mult_helper(v, self.left_multiplication_matrices)

    def right_multiplication_matrix(self,v):
        return self.mult_helper(v, self.right_multiplication_matrices)


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

        # useful optimization for implementing multiply
        self.multiplication_matrix = sparse.hstack(left_multiplication_matrices)

        # unit of algebra. Dense numpy array
        self.unit = unit


class Algebra:

    def associative_defect(self):
        m = self.multiplication
        return np.sum(np.abs(np.tensordot(m,m,(2,0)) -
                             np.transpose(np.tensordot(m,m,(1,2)),(0,2,3,1))))

    def left_identity_defect(self):
        m = self.multiplication
        return np.sum(
            np.abs(np.tensordot(self.unit, m, (0,0)) -
                   np.identity(self.dimension)))

    def right_identity_defect(self):
        m = self.multiplication
        return np.sum(
            np.abs(np.tensordot(self.unit, m, (0,1)) -
                   np.identity(self.dimension)))


    def commutative_defect(self):
        m = self.multiplication
        return np.sum(np.abs(m - np.transpose(m,(1,0,2))))

    def algebra_defect(self):
        return sum([self.associative_defect(),
                    self.left_identity_defect(),
                    self.right_identity_defect()])

    def irrep_defect_multiplication(self,irrep):
        m = self.multiplication
        result = np.sum(np.abs(np.tensordot(m, irrep, (2,0)) -
                               np.transpose(np.tensordot(irrep, irrep, (2,1)),(2,0,1,3))))
        return result

    def irrep_defect_identity(self,irrep):
        m = self.multiplication
        result = np.sum(np.abs(np.tensordot(self.unit, irrep, (0,0)) - np.identity(irrep.shape[1])))
        return result

    def irrep_defect(self,irrep):
        return self.irrep_defect_multiplication(irrep) + self.irrep_defect_identity(irrep)

    def multiply(self,x,y):
        m = self.multiplication
        return np.tensordot(np.tensordot(x,m,(0,0)),y,(0,0))

    def random_vector(self):
        d = self.dimension
        return (np.random.rand(d) - 0.5 ) + 1j * (np.random.rand(d) - 0.5 )

    def commutator_matrix(self):
        d = self.dimension
        m = self.multiplication
        return np.transpose(np.reshape(m - np.transpose(m,(1,0,2)),(d,d*d)))

    def left_multiplication_matrix(self,v):
        m = self.multiplication
        return np.transpose(np.tensordot(v, m, (0,0)))

    def right_multiplication_matrix(self,v):
        m = self.multiplication
        return np.transpose(np.tensordot(v,m,(0,1)))

    def __init__(
            self,
            dimension,
            multiplication,
            unit,
            is_sparse = False
    ):

        # dimension of algebra
        # we denote the basis vectors as e_1, e_2, ...., e_dimension
        self.dimension = dimension

        # multiplication tensor (3D numpy array)
        # e_i e_j = sum over k multiplication[i,j,k] e_k
        self.multiplication = multiplication

        # unit of the algebra in the given basis
        self.unit = unit

        self.is_sparse = is_sparse



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

