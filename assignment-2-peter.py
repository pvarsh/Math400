"""
Math 400
Assignment 2
Author: Peter Varshavsky
"""

import matplotlib.pyplot as plt
from math import *
from operator import mul

from scipy.interpolate import interp1d

"""
Problem 3
"""

"""
Function definitions
"""

def listProd(X):
    return reduce(mul, X)

def lagrange_i(X, x, i):
    index = range(len(X))
    index.remove(i)

    denom = listProd([X[i] - X[j] for j in index])
    numer = listProd([x - X[j] for j in index])
    
    return numer/denom

def lagrange(X, Y, x):

    N = len(X)
    return sum([Y[i]*lagrange_i(X, x, i) for i in range(N)])

def stupidSearch(X, x):
    "Returns the index of greatest element smaller than x"

    if x < X[0] or x > X[-1]:
        print "stupidSearch: search outside range of values"
        return -1
    for i in range(1, len(X)):
        if X[i] > x:
            return i-1
    return i
    


def makeCubicSpline(X, Y):
    N = len(X) - 1

    # Copy Y
    A = [y for y in Y]

    print "N: ", N
    print "len(X): ", len(X)
    
    # Step 1
    H = [X[i+1] - X[i] for i in range(N)]
    print "len(H): ", len(H)
    # Step 2
    Alpha = [3.0/H[i] * (A[i+1] - A[i]) - 3.0/H[i-1] * (A[i] - A[i-1]) for i in range(1, N)]
    print "len(Alpha): ", len(Alpha)
    
    # Step 3
    L = [1]
    Mu = [0]
    Z = [0]

    # Step 4
    for i in range(1, N):
        L.append(2 * (X[i+1] - X[i]) - H[i-1]*Mu[i-1])
        Mu.append(H[i] / L[i])
        Z.append( (Alpha[i-1] - H[i-1] * Z[i-1] ) / L[i] )

    # Step 5
    L.append(1)
    Mu.append(0)
    C = [0 for _ in range(N+1)]
    B = [0 for _ in range(N+1)]
    D = [0 for _ in range(N+1)]

    # Step 6
    for j in reversed(range(N)):
        print "j: ", j
        C[j] = Z[j] - Mu[j] * C[j+1]
        B[j] = (A[j+1] - A[j]) / H[j] - H[j] * (C[j+1] + 2* C[j]) / 3
        D[j] = (C[j+1] - C[j]) / (3*H[j])

    print A
    print B
    print C
    print D
    
    return A, B, C, D


def cubicSplineInterpolate(X, x, A, B, C, D):
    j = stupidSearch(X, x)
    if j == -1:
        print "x value outside data"
        return -1

    return A[j] + B[j]*(x - X[j]) + C[j] * (x - X[j])**2 + D[j] * (x - X[j])**3
        


##def lagrange_fun_maker(X, Y, x):
##    N = len(X)
##    def temp(x):
##        return sum([Y[i]*lagrange_i(X, x, i) for i in range(N)])
##    return temp

##def binary_search(X, x):
##    N = len(X)
##    if x == X[N/2]:
##        return


##def linearSpline(X, Y, x):
##    if x > X[-1] or x < X[0]:
##        print "Warning, this is an extrapolation"
##
##    for i, j in enumerate(X):
##        if j <= x


    
def f_prob3(x):
    return 1.6*exp(-2*x) * sin(3*pi*x)

X = [0, 1.0/6, 1.0/3, 1.0/2, 7.0/12, 2.0/3, 3.0/4, 5.0/6, 11.0/12, 1]
Y = [1.6 * exp(-2 * x) * sin(3*pi*x) for x in X]
Yprime = [-3.2 * exp(-2*x) * sin(3*pi*x) + 4.8 * exp(-2*x) * cos(3*pi*x) for x in X]



A, B, C, D = makeCubicSpline(X, Y)

#print Y
#print Yprime

how_fine_Orig = 100
how_fine_Lagr = 100
X_fine_Orig = [X[0] + (float(X[-1]) - X[0])/how_fine_Orig * i for i in range(how_fine_Orig + 1)]
X_fine_Lagr = [X[0] + (float(X[-1]) - X[0])/how_fine_Lagr * i for i in range(how_fine_Lagr + 1)] 
Y_fine_Lagr = [lagrange(X, Y, x) for x in X_fine_Lagr]
Y_fine_Orig = [f_prob3(x) for x in X_fine_Orig]
Y_fine_CS = [cubicSplineInterpolate(X, x, A, B, C, D) for x in X_fine_Orig]

f2 = interp1d(X, Y, kind='cubic')

plt.plot(X, Y, "rs")
plt.plot(X_fine_Orig, f2(X_fine_Orig), color = "red", label = "Scipy cubic")
plt.plot(X_fine_Lagr, Y_fine_Lagr, color = "orange", label = "Lagrange")
plt.plot(X_fine_Orig, Y_fine_Orig, color = "cyan", label = r"$1.6e^{-2x}\sin(3\pi x)$")
plt.plot(X_fine_Orig, Y_fine_CS, color = "magenta", label = "Cubic Spline")
plt.legend(loc = "lower right")
plt.show()





