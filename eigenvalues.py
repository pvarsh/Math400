"""
Finding eigenvalues by power mehtod
Author: Peter Varshavsky
"""

from gaussian_elimination_2 import *
from random import randint

def infty_norm(vec):
    return max(vec)[0]


def eigen_power(A, tolerance, maxIter):
    x0 = [randint(2, 10) for i in range(rows(A))]

    for i in range(maxIter):
        x = matMult(A, transpose([x0]))
        print "norm"
        print infty_norm(x)
        x = scalarMult(1.0/infty_norm(x), x)

        print "x: ", x
        print "x0: ", x0

        if infty_norm(subVectors(x, x0)) < tolerance:
            print("found")
            return x

        x0 = x

A = [[3, 1, 0], [1, 2, 3], [1, 2, 10]]
eigen_power(A, 0.1, 3)
