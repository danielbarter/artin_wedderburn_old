# ArtinWedderburn

A package for numerically computing irreducible representations of semi-simple algebras over the complex numbers.

### description

Take a semi-simple algebra with basis `e1, ..., ed`. The multiplication tensor is defined by
```
ei ej = sum over k m[i,j,k] ek
```
The algebra unit vector is defined by
```
ei u = u ei = ei for each i
```
ArtinWedderburn takes a multiplication tensor and unit vector and computes the irreducible representations for the algebra (these are only defined upto conjugation). It keeps track of numerical error as it works, so you can see if it is polluting the result.

### usage

You can compute the irreducible representations of S4 as follows:

```
from ArtinWedderburn import *
alg = load_algebra_from_file('./symmetric_groups/4')
aw = ArtinWedderburn(alg, logging=True)
```

if you are looking for something more exotic, Jacob Bridgeman has computed the tube algebras for all multiplicity free unitary fusion categories of rank < 7: [10.5281/zenodo.4277499](https://zenodo.org/record/4277499). There are 195 such algebras in total with dimensions ranging from 4 to 301. The algebras are stored in sparse format in the folder `small_tube_algebras`. All of their irreducible representations can be computed like this:

```
python benchmark.py small_tube_algebras
```

This runs in 2 minutes with a max ram usage of 328MB on a Intel Xeon W-1290 CPU. We have also compute the irreps of some larger algebras on a Intel Xeon E5-2670 CPU:

| dimension  | time (mins) | max ram usage (GB) |
|     ---    |     ----    |       ---          |
| 783 | 16 | 4 |
| 910 | 22 | 7 |
| 1112 | 93 | 10 |
| 1547 | 517 | 24 |

Each of these algebras has 6 irreps.

### dependencies

ArtinWedderburn uses `numpy` for tensor arithmetic and `scipy` for sparse matrices, SVDs and computing eigenvalues. It should work with any recent versions of those libraries.
