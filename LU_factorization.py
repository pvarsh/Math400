"""
LU factorization
Author: Peter Varshavsky
LUfactorize() is the exact implementation of LU Factorization algorithm from
    page 406 of Numerical Analysis 9th ED by Richard Burden
Started: 3/10/14
"""

from gaussian_elimination_2 import *
from lu import *
from mirai_LU import lUFactorization

def partialDotProduct(vec1, vec2, kmax):
    return dot(vec1[0:kmax+1], vec2[0:kmax+1])

def LUfactorize_(A):
    
    n = rows(A) - 1
    L = identity(n + 1)
    U = zero(n+1,n+1)

    # Step 1
    print("Step 1")
    U[0][0] = A[0][0]
    if A[0][0] == 0:
        print "Factorization impossible"
        return None

    # Step 2
    print("Step 2")
    for j in range(1, n+1):
        U[0][j] = float(A[0][j])/L[0][0]
        print "\nU: "
        show(U)
        L[j][0] = float(A[j][0])/U[0][0]
        print "\nL: "
        show(L)

    # Step 3
    print("Step 3")
    for i in range(1, n):
        # Step 4
        print("Step 4")
        U[i][i] = float((A[i][i] - partialDotProduct(getRow(L, i), getCol(U, i), i-1)))/L[i][i]
##      print "Weird error in LUfactorize_"
##      print "Factorization impossible or something"
##      return(-1)
        if U[i][i]*L[i][i] == 0:
            print "Factorization impossible"
            return None
        else:
            # Step 5
            print("Step 5")
            for j in range(i+1, n+1):
                U[i][j] = (1.0 / L[i][i]) * (A[i][j] - partialDotProduct(getRow(L, i), getCol(U, j), i-1)) 
                L[j][i] = (1.0 / U[i][i]) * (A[j][i] - partialDotProduct(getRow(L, j), getCol(U, i), i-1))
    # Step 6
    print("Step 6")
    U[n][n] = A[n][n] - partialDotProduct(getRow(L, n), getCol(U, n), n)

    # Step 7
    return [L, U]

def solve_LU(A, b):

    [L, U] = LUfactorize(A)
    show(L)

    y = [ float(b[0])/L[0][0] ]

    for i in range(1, rows(A)):
        y.append( (1.0 / L[0][0]) * (b[i] - partialDotProduct(L[i], y, i-1)))

    x = backSub(augment(U, y))
    
    print "x is: "
    print x
    return x


    


A = [[1,1,0,3],
     [2,1,-1,1],
     [3,-1,-1,2],
     [-1,2,3,-1]]

b = [8,7,14,-7]

#solve_LU(A, b)

def test_LUfactorize():
        
    A = [[1,1,0,3],
         [2,1,-1,1],
         [3,-1,-1,2],
         [-1,2,3,-1]]
    L, U = LUfactorize_(A)
    print "Finished factorization"
    print "A: "
    show(A)
    print "L: "
    show(L)
    print "U: "
    show(U)
    print "L * U: "
    show(matMult(L, U))

#test_LUfactorize()


def testLU():

    from random import randint
    from random import seed

    seed(1234)
    n = 4
    A = [[randint(1, 20) for j in range(n)] for i in range(n)]
    b = [randint(1, 20) for i in range(n)]
    A = [[1,1,0,3],
     [2,1,-1,1],
     [3,-1,-1,2],
     [-1,2,3,-1]]
    print "A"
    show(A)
    print "b"
    print(b)

    print "My LU factorization"
    [L1, U1] = LUfactorize(A)
    A_ = matMult(L1, U1)
    show(A_)

    print "Lu factorization from internet (lu.py)"
    [_, L2, U2] = lu_decomposition(A)
    A_ = mult_matrix(L2, U2)
    show(A_)

    print "Mirai LU"
    [L3, U3] = lUFactorization(A)
    A_ = matMult(L3, U3)
    show(A_)

    print "SciPy LU"
    import pprint
    import scipy
    import scipy.linalg   # SciPy Linear Algebra Library

    A = scipy.array(A)
    P, L4, U4 = scipy.linalg.lu(A)
    A_ = matMult(L4, U4)
    show(A_)
    
    
    
    print "My L: "
    show(L1)
    print "Internet L: "
    show(L2)
    print "Mirai L: "
    show(L3)
    print "My U: "
    show(U1)
    print "Internet U: "
    show(U2)
    print "Mirai U: "
    show(U3)
    

    
    """
    print "************************\n************************"

    solLU = solve_LU(A, b)
    solGE = ge_1(augment(A, b))[1]

    print "solution LU: ", solLU
    print "solution GE: ", solGE

    print "residueTest LU: "
    residueTest(A, b, solLU)
    print "residueTest GE: "
    residueTest(A, b, solGE)
    """

testLU()
"""
A = [[1,2,    3],
     [4,   5,    6],
     [7,   8,    0 ]]
[L, U] = LUfactorize(A)
show(L)
show(U)
show(matMult(L, U))
"""
