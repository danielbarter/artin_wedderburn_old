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
ArtinWedderburn takes a multiplication tensor and unit vector and computes the irreducible representations for the algebra (these are only defined upto conjugation). It keeps track of numerical error as it works, so you can see if it is polluting the result. Currently, it can compute the irreducible representations for S5 in < 1 minute and 3GB of memory on a Ryzen 5 2400G.

### usage

```
from ArtinWedderburn import *
alg = symmetric_group_algebra(3)
aw = ArtinWedderburn(alg, logging=True)
```

Algebra objects can be constructed using `Algebra(dimension, multiplication, unit)` where `multiplication` is a numpy array with shape `(dimension,dimension,dimension)` and `unit` is a numpy vector with length `dimension`.

### dependencies

ArtinWedderburn uses `numpy` for tensor arithmetic and `scipy` for SVDs and computing eigenvalues. It should work with any recent versions of those libraries.
