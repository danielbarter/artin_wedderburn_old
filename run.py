import sys
import os
from ArtinWedderburn import *
import cProfile, pstats, io
from pstats import SortKey


if len(sys.argv) != 2:
    print("usage: python run.py algebra_file")
    exit()

with cProfile.Profile() as pr:
    alg = load_algebra_from_file(sys.argv[1])
    aw = ArtinWedderburn(alg, logging=True)


s = io.StringIO()
sortby = SortKey.CUMULATIVE
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()
print(s.getvalue())
