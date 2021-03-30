import sys
import os
from ArtinWedderburn2 import *


if len(sys.argv) != 2:
    print("usage: python benchmark.py algebra_folder")
    exit()

directory = sys.argv[1]
os.chdir(directory)
algebra_files = os.listdir()

for path in algebra_files:
    alg = load_sparse_algebra_from_file(path)
    aw = ArtinWedderburn(alg)
    print(path)
    print("algebra dimension = ", alg.dimension)

    print("irrep dimensions =", end = ' ')
    for dim in aw.irrep_dimensions.values():
            print(dim, end=' ')


    print("\n",end='')
    print("total defect = ", aw.total_defect)
    print("\n")

