"""
Hilbert Matrix Exercise
Code: Peter Varshavsky
Started: 03/10/14
"""

from gaussian_elimination_2 import show

def makeHilbert(n):
    " Returns an nxn Hilbert matrix "
    " Author: PV "
    return [[1.0/(col+1) for col in range(row, n+row)] for row in range(0, n)]

def integralRecursion(n):
    " Returns a vector of $\int_0^1 x^k \sin(\pi x)$ for k from 0 to n "
    " Author: PV "
    from math import pi
    b = [2/pi]
    b.append(1/pi)

    for i in range(2, n):
        b.append(1/pi - i*(i-1)*b[-2]/(pi**2))
    return b


