import sys
sys.path.append("../")

from ArtinWedderburn import *

with open("test.txt",'r') as f:
    lines = f.readlines()

identity_vec = [int(x) for x in lines[0].split(' ')]
dim = len(identity_vec)
unit_coefficients = []
for (index, val) in enumerate(identity_vec):
    if val != 0:
        unit_coefficients.append(((index,),float(val)))

unit = sparse.COO.from_iter(unit_coefficients, shape = (dim,), dtype=complex)

multiplication_lines = lines[1:]
multiplication_coefficients = []
for line in multiplication_lines:
    parts = line.split(' ')
    i = int(parts[0]) - 1
    j = int(parts[1]) - 1
    k = int(parts[2]) - 1
    val = float(parts[3])
    multiplication_coefficients.append(((i,j,k),val ))

mult = sparse.COO.from_iter(multiplication_coefficients, shape=(dim,dim,dim), dtype=complex)

alg = Algebra(dim, mult, unit, is_sparse=True)
