import numpy as np
import scipy.sparse as sparse
from scipy.linalg import eigh, eig, svd
from math import factorial, sqrt

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


def eigenspace(matrix):
    eigvals, eigvecs = eigh(np.dot(matrix, np.transpose(matrix.conjugate())))
    u = np.flip(eigvecs, 1)
    s = np.flip(eigvals, 0)
    return u, s


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
            algebra.multiply(center_inclusion[:,i], center_inclusion[:,j]))
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
        idempotent = self.unscaled_central_idempotents[:, idempotent_index]
        left_multiplication = algebra.left_multiplication_matrix(idempotent)

        # u is the change of basis from new basis to old basis
        # u, s, v = svd(left_multiplication)
        u, s = eigenspace(left_multiplication)

        # u is unitary, so inverse is conjugate transpose
        # u_inverse is change of basis from old basis to new basis
        # u_inverse = np.transpose(np.conj(u))

        # the basis for the block is the columns of u whose singular value
        # is above the threshold

        counter = 0
        for i in range(len(s)):
            if abs(s[i]) > self.threshold:
                counter += 1

        block_inclusion = u[:,:counter]
        self.block_inclusions[idempotent_index] = block_inclusion
        block_projection = np.conj(np.transpose(block_inclusion))
        block_dimension = block_inclusion.shape[1]


        block_multiplication = np.array([[ np.dot(
            block_projection,
            algebra.multiply(block_inclusion[:,i], block_inclusion[:,j]))
           for i in range(block_dimension)]
         for j in range(block_dimension)])

        block_pre_unit = np.dot(block_projection, idempotent)

        pre_block = Algebra(
            block_dimension,
            block_multiplication,
            block_pre_unit)

        # we have only captured the identity up to scalar multiplication
        # we need the identity on the node if we want to recover the correct scale
        # for the representations
        factor = pre_block.left_multiplication_matrix(pre_block.unit)[0,0]
        block = Algebra(
            pre_block.dimension,
            pre_block.multiplication,
            pre_block.unit / factor)

        block_defect = block.algebra_defect()
        self.log("block defect:", format_error(block_defect))
        self.total_defect += block_defect

        central_idempotent = np.dot(block_inclusion,block.unit)


        self.blocks[idempotent_index] = block
        self.central_idempotents[idempotent_index] = central_idempotent


    def compute_irrep(self,block_index):

        block = self.blocks[block_index]
        sqrt_dimension = int(sqrt(block.dimension))
        self.irrep_dimensions[block_index] = sqrt_dimension


        inclusion = self.block_inclusions[block_index]
        projection = np.transpose(np.conj(inclusion))

        if block.dimension == 1:
            block_irrep = block.multiplication
            m_defect = 0.0

        else:

            v = block.random_vector()
            vr = block.right_multiplication_matrix(v)
            eigenvalues_with_repeats = eig(vr)[0]
            eigenvalues = fuzzy_filter(eigenvalues_with_repeats,self.threshold)

            accumulator = block.unit
            for i in range(sqrt_dimension - 1):
                accumulator = block.multiply(accumulator, v - eigenvalues[i] * block.unit)

            proj = block.right_multiplication_matrix(accumulator)

            # u goes from new basis to old basis
            # u, s, v = svd(proj)
            u, s = eigenspace(proj)

            # u_inverse goes from old basis to new basis
            u_inverse = np.transpose(np.conj(u))


            m = block.multiplication
            change_codomain_basis = np.tensordot(m,u_inverse, (2,1))
            m_rotated = np.transpose(np.tensordot(u,change_codomain_basis,(0,1)),(1,0,2))
            m_defect = np.sum(np.abs(m_rotated[:,
                                               0:sqrt_dimension,
                                               sqrt_dimension:block.dimension]))

            block_irrep = m_rotated[:,0:sqrt_dimension,0:sqrt_dimension]

        self.log("SVD defect:", format_error(m_defect))
        self.total_defect += m_defect

        block_irrep_defect = block.irrep_defect(block_irrep)
        self.log("block irred defect:", format_error(block_irrep_defect))
        self.total_defect += block_irrep_defect

        irrep = np.tensordot(projection, block_irrep, (0,0))
        self.irreps[block_index] = irrep

        irrep_defect = self.algebra.irrep_defect(irrep)
        self.log("irrep defect:", format_error(irrep_defect))
        self.total_defect += irrep_defect

    def central_idempotent_defect(self):
        center = self.center
        algebra = self.algebra
        central_idempotents = self.central_idempotents
        defect_mat = np.zeros((center.dimension, center.dimension))
        for i in range(center.dimension):
            for j in range(center.dimension):
                if i != j:
                    defect_mat[i,j] = np.sum(np.abs(algebra.multiply(
                        central_idempotents[i],
                        central_idempotents[j])))
                else:
                    defect_mat[i,j] = np.sum(np.abs(algebra.multiply(
                        central_idempotents[i],
                        central_idempotents[j]) - central_idempotents[i]))

        id_defect = np.copy(algebra.unit)

        for i in range(center.dimension):
            id_defect -= central_idempotents[i]

        return np.sum(defect_mat) + np.sum(np.abs(id_defect))


    def central_idempotent_irrep_defect(self):
        center = self.center
        central_idempotents = self.central_idempotents
        irreps = self.irreps

        defect_mat = np.zeros((center.dimension, center.dimension))

        for i in range(center.dimension):
            for j in range(center.dimension):
                idempotent = central_idempotents[i]
                irrep = irreps[j]
                if i != j:
                    defect_mat[i,j] = np.sum(np.abs(np.tensordot(
                        idempotent,
                        irrep,
                        (0,0))))

                else:
                    mat = np.tensordot(
                        central_idempotents[i],
                        irreps[j],
                        (0,0))

                    mat -= np.identity(irrep.shape[-1])

                    defect_mat[i,j] = np.sum(np.abs(mat))

        return np.sum(defect_mat)



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

        self.irreps = {}
        self.irrep_dimensions = {}
        self.log("computing irreps...")
        self.log("")

        for i in range(self.center.dimension):
            self.log("computing irrep", i)
            self.compute_irrep(i)
            self.log("")

        self.log("all done!")

        central_idempotent_defect = self.central_idempotent_defect()
        self.log("central idempotent defect:", format_error(central_idempotent_defect))
        self.total_defect += central_idempotent_defect

        central_idempotent_irrep_cross_defect = self.central_idempotent_irrep_defect()
        self.log("central idempotent irrep cross defect:",
                 format_error(central_idempotent_irrep_cross_defect))
        self.total_defect += central_idempotent_irrep_cross_defect


        self.log("total defect:", format_error(self.total_defect))

        self.log("irrep tensors are stored in the attribute self.irreps")
        self.log(self.center.dimension, "irreducible representations with dimensions")
        for index, dim in self.irrep_dimensions.items():
            self.log(index, ":", dim)




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
        xy = self.multiply(x,y)

        # irrep[i,j] is the vector e_i b_j
        l = np.tensordot(xy, irrep, (0,0))

        r_uneval = np.tensordot(irrep, irrep, (2,1))

        # this is backwards because irrep[i] is the transpose of
        # what would conventionally be called the action matrix
        r = np.tensordot(y, np.tensordot(x, r_uneval, (0,0)), (0,1))
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
        return self.multiplication_matrix * np.kron(x,y)

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
        x = self.random_vector()
        y = self.random_vector()
        z = self.random_vector()
        l = self.multiply(self.multiply(x,y), z)
        r = self.multiply(x, self.multiply(y,z))
        return np.sum(np.abs(l - r))



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
                    self.associative_defect(),
                    self.associative_defect(),
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
            unit
    ):

        # dimension of algebra
        # we denote the basis vectors as e_1, e_2, ...., e_dimension
        self.dimension = dimension

        # multiplication tensor (3D numpy array)
        # e_i e_j = sum over k multiplication[i,j,k] e_k
        self.multiplication = multiplication

        # unit of the algebra in the given basis
        self.unit = unit




def load_algebra_from_file(path):
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

