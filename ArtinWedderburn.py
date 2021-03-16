import numpy as np

from itertools import permutations
from math import factorial, sqrt
from scipy.linalg import svd,eig

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
        m = algebra.multiplication
        unit = algebra.unit
        d = algebra.dimension


        u,s,v = svd(algebra.commutator_matrix())
        # v is the change of basis matrix from our old basis to our new basis
        # v is unitary, so inverse is conjugate transpose
        # v_inverse is the change of basis from our new basis to our old basis
        v_inverse = np.transpose(np.conj(v))

        # the singular values in s are in non increasing order
        # this means that the center basis vectors will the last columns of v_inverse
        # whrere singular values are below threshold

        counter = 0
        for i in range(len(s)):
            if abs(s[i]) < self.threshold:
                counter += 1

        center_inclusion = v_inverse[:,d - counter:]
        self.center_inclusion = center_inclusion
        center_dimension = center_inclusion.shape[1]


        change_codomain_basis = np.tensordot(m, v, (2,1))
        change_left_basis = np.tensordot(v_inverse,
                                         change_codomain_basis,
                                         (0,0))

        m_rotated_anti = np.tensordot(
            v_inverse,
            change_left_basis,
            (0,1))

        m_rotated = np.transpose(m_rotated_anti,(1,0,2))

        m_defect = np.sum(np.abs(m_rotated[d-center_dimension:d,
                                                d-center_dimension:d,
                                                0:d-center_dimension]))

        self.log("SVD multiplication defect:",format_error(m_defect))
        self.total_defect += m_defect

        unit_rotated = np.tensordot(unit, v, (0,1))

        unit_defect = np.sum(np.abs(unit_rotated[0:d-center_dimension]))
        self.log("SVD unit defect:", format_error(unit_defect))
        self.total_defect += unit_defect

        center = Algebra(center_dimension,
                         m_rotated[d-center_dimension:d,
                                   d-center_dimension:d,
                                   d-center_dimension:d],
                         unit_rotated[d-center_dimension:d])

        self.center = center

        center_defect = center.algebra_defect()
        self.log("center defect:", format_error(center_defect))
        self.total_defect += center_defect


        center_commutative_defect = center.commutative_defect()
        self.log("center commutative defect:", format_error(center_commutative_defect))
        self.total_defect += center_commutative_defect


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


        pre_block = Algebra(
            block_dimension,
            m_rotated[0:block_dimension,
                      0:block_dimension,
                      0:block_dimension],
            unit_rotated[0:block_dimension])

        # we have only captured the identity upto scalar multiplication
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
            u, s, v = svd(proj)

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


        irrep = np.tensordot(projection, block_irrep, (0,0))
        self.irreps[block_index] = irrep

        irrep_defect = self.algebra.irrep_defect(irrep)
        self.log("irrep defect:", format_error(irrep_defect))
        self.total_defect += irrep_defect


    #TODO: write a function to validate the irreps using the central idempotents

    def log(self,*args, **kwargs):
        if self.logging:
            print(*args, **kwargs)

    def __init__(self, algebra, threshold = 1.0e-5, logging = False):
        self.algebra = algebra
        self.threshold = threshold
        self.logging = logging
        self.total_defect = algebra.algebra_defect()

        self.log("algebra defect:", self.total_defect)
        self.log("")

        self.log("computing center...")
        self.compute_center()
        self.log("")

        # really, this only computes the central idempotents upto scalar multiplication
        # we compute the central idempotents on the nose while computing the blocks
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
        self.log("total defect:", format_error(self.total_defect))

        self.log("irrep tensors are stored in the attribute self.irrep")
        self.log(self.center.dimension, "irreducible representations with dimensions")
        for index, dim in self.irrep_dimensions.items():
            self.log(index, ":", dim)



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
            unit):

        # dimension of algebra
        # we denote the basis vectors as e_1, e_2, ...., e_dimension
        self.dimension = dimension

        # multiplication tensor (3D numpy array)
        # e_i e_j = sum over k multiplication[i,j,k] e_k
        self.multiplication = multiplication

        # unit of the algebra in the given basis
        self.unit = unit


# we represent a permutation as an n-tuple containing the numbers 0,1,...,n-1
# if we think of such a tuple as a function from indices to values
# then permutations are composed via function composition as follows
def multiply_permutations(p2, p1):
    return tuple([p2[i] for i in p1])

def symmetric_group_algebra(n):
    dimension = factorial(n)
    multiplication = np.zeros((dimension,dimension,dimension), dtype=complex)
    unit = np.zeros(dimension,dtype=complex)


    permutation_to_index = {}
    index_to_permutation = {}
    index = 0



    for permutation in permutations(range(n)):
        permutation_to_index[permutation] = index
        index_to_permutation[index] = permutation
        index += 1


    identity_permutation = tuple(range(n))
    unit[permutation_to_index[identity_permutation]] = 1.0

    for i in range(dimension):
        for j in range(dimension):
            pi = index_to_permutation[i]
            pj = index_to_permutation[j]
            pk = multiply_permutations(pi,pj)
            multiplication[i,j,permutation_to_index[pk]] = 1.0

    return Algebra(dimension, multiplication, unit)
