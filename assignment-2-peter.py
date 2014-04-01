"""
Math 400
Assignment 2
Author: Peter Varshavsky
"""

import matplotlib.pyplot as plt
from math import *
from operator import mul
from scipy.interpolate import interp1d # using this to compare with my interpolations

"""
Problem 3
"""

"""
Function definitions
"""


def listProd(X):
    "Returns the product of items in a list"
    return reduce(mul, X)

def lagrange_i(X, x, i):
    "Returns L_i(x) using the product form formula from exercise 1"
    index = range(len(X))
    index.remove(i) # this makes sure denominator and numerator don't have (x_i - x_i)

    denom = listProd([X[i] - X[j] for j in index])
    numer = listProd([x - X[j] for j in index])
    
    return numer/denom

def lagrange(X, Y, x):
    "Returns the Lagrange interpolation value at x"
    "Uses lagrange_i()"
    N = len(X)
    return sum([Y[i]*lagrange_i(X, x, i) for i in range(N)])

def stupidSearch(X, x):
    "Returns the index of greatest element smaller than x"
    "Called stupid because it searches linearly, which is inefficient"

    if x < X[0] or x > X[-1]:
        print "stupidSearch: search outside range of values"
        return -1
    for i in range(1, len(X)):
        if X[i] > x:
            return i-1
    return i
    


def makeCubicSpline(X, Y):
    "Implementation of natural cubic spline algorithm from textbook"
    "Returns four lists of coefficients A, B, C, D"
    N = len(X)

    # Copy Y
    A = [y for y in Y]

    #print "N: ", N
    #print "len(X): ", len(X)
    
    # Step 1
    H = [ X[i+1] - X[i] for i in range(N-1) ]
    #print "len(H): ", len(H)

    # Step 2
    Alpha = [3.0/H[i] * (A[i+1] - A[i]) - 3.0/H[i-1] * (A[i] - A[i-1]) for i in range(1, N-1)]
    #print "len(Alpha): ", len(Alpha)
    
    # Step 3
    L = [1]
    Mu = [0]
    Z = [0]

    # Step 4
    for i in range(1, N-1):
        L.append(2 * (X[i+1] - X[i-1]) - H[i-1]*Mu[i-1])
        Mu.append(H[i] / L[i])
        Z.append( (Alpha[i-1] - H[i-1] * Z[i-1] ) / L[i] ) # Alpha is i-1 since in the book
                                                           # it is defined with 1-based index

    # Step 5
    L.append(1)
    Z.append(0)
    C = [0 for _ in range(N+1)]
    B = [0 for _ in range(N+1)]
    D = [0 for _ in range(N+1)]

    # Step 6
    for j in reversed(range(N-1)):
        #print "j: ", j
        C[j] = Z[j] - Mu[j] * C[j+1]
        B[j] = (A[j+1] - A[j]) / H[j] - H[j] * (C[j+1] + 2 * C[j]) / 3
        D[j] = (C[j+1] - C[j]) / (3*H[j])

    #print A
    #print B
    #print C
    #print D
    
    return A, B, C, D


def cubicSplineInterpolate(X, x, A, B, C, D):
    "Returns the cubic spline interpolation value at x"
    "Coefficients A, B, C, D need to be computed by makeCubicSpline()"
    j = stupidSearch(X, x)
    if j == -1:
        print "x value outside data"
        return -1

    return A[j] + B[j]*(x - X[j]) + C[j] * (x - X[j])**2 + D[j] * (x - X[j])**3
        
def linearSplineInterpolate(X, Y, x):
    "Returns the interpolation value at x using linear splines"
    
    j = stupidSearch(X, x)

    if j == -1:
        print "x value outside data"
        return -1

    if j == len(X) - 1:
        return Y[-1]
    
    return Y[j] - Y[j] * (x - X[j]) / (X[j+1] - X[j]) + Y[j+1] * (x - X[j]) / (X[j+1] - X[j])

def l_i_prime(X, i):
    index = range(len(X))
    index.remove(i)
    print "\nindex in l_i_prime: ", index

    return l_i_prime_recursion(X, i, 0, index)
    

def l_i_prime_recursion(X, i, j, index):
    if j == len(index) - 2:
        return (X[i] - X[index[j]]) + X[i] - X[index[-1]]
    else:
        first_term = X[i] - X[index[j]]
##def l_i_prime_recursion(X, i, j, index):
##
##    if j == len(index)-1:
##        print "reached recursion end"
##        return (X[i] - X[-1]) + (X[i] - X[j])
##    else:
##        print "recursion:: j = ", j
##        print "index: ", index[j: ]
##        
##        first_summand_list = [ X[i] - X[k] for k in index[j: ] ]
##        print "summing over index: ", index[j+1: ]
##        print "first summand list: ", first_summand_list
##        first_summand = listProd(first_summand_list)
##        print "first summand: ", first_summand
##        print "X[i] - X[j]: ", X[i], X[j]
##
##        return first_summand + (X[i]-X[j]) * l_i_prime_recursion(X, i, j+1, index)
##                    
    
def hermite_1_j(X, x, j):
    return (1 - 2*(x - X[j]) * l_i_prime(X, j) ) * lagrange_i(X, x, j)**2

    
def hermite_2_j(X, x, j):
    return (x - X[j]) * lagrange_i(X, x, j)**2

def hermite(X, Y, Yprime, x):
    return sum( [Y[j]*hermite_1_j(X, x, j) + Yprime[j]*hermite_2_j(X, x, j) for j in range(len(X))])

##def lagrange_fun_maker(X, Y, x):
##    N = len(X)
##    def temp(x):
##        return sum([Y[i]*lagrange_i(X, x, i) for i in range(N)])
##    return temp


    
def f_prob3(x):
    "Returns f(x), f(x) = 1.6 & exp(-2*x) * sin(3*pi*x)"
    return 1.6*exp(-2*x) * sin(3*pi*x)

"""
End function definitions
"""

"""
Function calls
"""


def exercise3():
    X = [0, 1.0/6, 1.0/3, 1.0/2, 7.0/12, 2.0/3, 3.0/4, 5.0/6, 11.0/12, 1]
    Y = [1.6 * exp(-2 * x) * sin(3*pi*x) for x in X]
    Yprime = [-3.2 * exp(-2*x) * sin(3*pi*x) + 4.8 * exp(-2*x) * cos(3*pi*x) for x in X] # derivative for Hermitian

    A, B, C, D = makeCubicSpline(X, Y)

    # you can change this value to get a finer or coarser plot
    M = 200 # this is the number of points used to construct interpolation graphs

    # x values for interpolation
    X_inter = [X[0] + (float(X[-1]) - X[0])/M * i for i in range(M + 1)]

    Y_inter_H = [hermite(X, Y, Yprime, x) for x in X_inter]

    # y values with Lagrange
    Y_inter_L = [lagrange(X, Y, x) for x in X_inter]

    # y values with cubic splines
    Y_inter_CS = [cubicSplineInterpolate(X, x, A, B, C, D) for x in X_inter]

    # y values with linear splines
    Y_inter_LS = [linearSplineInterpolate(X, Y, x) for x in X_inter]

    # y values of the original function
    Y_f = [f_prob3(x) for x in X_inter]

    # y values using scipy cubic interpolation
    f2 = interp1d(X, Y, kind='cubic')


    plt.plot(X, Y, "rs")
    #plt.plot(X_inter, f2(X_inter), color = "red", label = "Scipy cubic")
    #plt.plot(X_inter, Y_inter_L, color = "orange", label = "Lagrange")
    plt.plot(X_inter, Y_f, color = "cyan", label = r"$1.6e^{-2x}\sin(3\pi x)$")
    #plt.plot(X_inter, Y_inter_CS, color = "magenta", label = "Cubic Spline")
    #plt.plot(X_inter, Y_inter_LS, color = "blue", label = "Linear spline")
    plt.plot(X_inter, Y_inter_H, color = "yellow", label = "Hermite")
    plt.legend(loc = "lower right")
    plt.show()

# exercise3()


def example_1():
    X = [1.3, 1.6, 1.9]
    Y = [0.6200860, 0.4554022, 0.2818186]
    Yprime = [-0.5220232, -0.5698959, -0.5811571]
    x = 1.5
    l_i = [50.0/9 * x**2 - 175.0/9 * x + 152.0/9,
           -100.0/9 * x**2 + 320.0/9 * x - 247.0/9,
           50.0/9 * x**2 - 145.0/9 * x + 104.0/9]
    l_i_prime_list = [100.0/9 * x - 175.0/9,
                 -200.0/9 * x + 320.0/9,
                 100.0/9 * x - 145.0/9]
    for i in range(len(X)):
        print "\nlagrange_i:    ", lagrange_i(X, x, i)
        print "L_i algebraic: ", l_i[i]
        print "l_i_prime:           ", l_i_prime(X, i)
        print "l_i_prime algebraic: ", l_i_prime_list[i]

#example_1()

def testing_recursion_quadratic():
    X = [1.2,1.6]
    derivatives = [2 * x  + (-X[0] - X[1]) for x in X] 
    print "derivatives: ", derivatives
    for i in range(len(X)):
        print l_i_prime(X, i)

#testing_recursion_quadratic()

def testing_recursion_cubic():
    X = [1.2, 1.5, 1.7]
    derivatives = [3 * x**2 + 2 * (-X[0]-X[1]-X[2])*x + (X[0]*X[1] + X[0]*X[2] + X[1]*X[2]) for x in X]
    print "derivatives: ", derivatives
    for i in range(len(X)):
        der = l_i_prime(X, i)
        print "%s'th derivative:           " %i, der
        print "%s'th derivative should be: " %i, derivatives[i]

# testing_recursion_cubic()

def testing_recursion_4():
    X = [1.1, 1.2, 1.5, 1.7]
    a,b,c,d = X
    x = 1.7
    print 4*x**3
    print 3*(a+b+c+d)*(x**2)
    print 2*(a*b+a*c+a*d+b*c+b*d+c*d)*x
    a*b*c + a*c*d + a*b*d + b*c*d
    
    derivatives = [4*x**3 - 3*(a+b+c+d)*x**2 + 2*(a*b+a*c+a*d+b*c+b*d+c*d)*x - (a*b*c + a*c*d + a*b*d + b*c*d) for x in X]
    print "derivatives: ", derivatives
    for i in range(len(X)):
        der = l_i_prime(X, i)
        print "%s'th derivative:           " %i, der
        print "%s'th derivative should be: " %i, derivatives[i]

testing_recursion_4()
