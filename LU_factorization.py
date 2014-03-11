"""
LU factorization
Author: Peter Varshavsky
Started: 3/10/14
"""

from gaussian_elimination_2 import *

def partialDotProduct(vec1, vec2, kmax):
    return dot(vec1[0:kmax], vec2[0:kmax])

def LUfactorize_(A):
    n = rows(A)
    L = identity(n)
    U = zero(n,n)

    # Step 1
    print("Step 1")
    U[0][0] = A[0][0]
    if A[0][0] == 0:
        print "Factorization impossible"
        return None

    # Step 2
    print("Step 2")
    for j in range(1, n):
        U[0][j] = float(A[0][j])/L[0][0]
        print "\nU: "
        show(U)
        L[j][0] = float(A[j][0])/U[0][0]
        print "\nL: "
        show(L)

    # Step 3
    print("Step 3")
    for i in range(1, n-1):
        

        # Step 4
        print("Step 4")

        U[i][i] = (A[i][i] - partialDotProduct(getRow(L, i), getCol(U, i), i-1))/L[i][i]
##      print "Weird error in LUfactorize_"
##      print "Factorization impossible or something"
##      return(-1)
        if U[i][i]*L[i][i] == 0:
            print "Factorization impossible"
            return(-1)
        else:
            # Step 5
            print("Step 5")
            for j in range(i+1, n):
                U[i][j] = (1.0 / L[i][i]) * (A[i][j] - partialDotProduct(getRow(L, i), getCol(U, j), i-1)) 
                L[j][i] = (1.0 / U[i][i]) * (A[j][i] - partialDotProduct(getRow(L, j), getCol(U, i), i-1)
)
    # Step 6
    print("Step 6 not written")

    # Step 7
    return([L, U])

#A = [[1,2,3],[1,2,4],[10,20,1]]
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
