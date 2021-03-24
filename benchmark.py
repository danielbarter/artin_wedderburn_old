import sys

if len(sys.argv) != 2:
    print("usage: python benchmark.py n")
    exit()



from ArtinWedderburn import *
alg = sparse_symmetric_group_algebra(int(sys.argv[1]))
aw = ArtinWedderburn(alg, logging=True)
