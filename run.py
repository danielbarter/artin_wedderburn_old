import sys
import os
from ArtinWedderburn import *


if len(sys.argv) != 2:
    print("usage: python run.py algebra_file")
    exit()


from ArtinWedderburn import *
alg = load_algebra_from_file(sys.argv[1])
aw = ArtinWedderburn(alg, logging=True)
