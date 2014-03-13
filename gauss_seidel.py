"""
Gauss - Seidel iterative method
Author: Peter Varshavsky
GS() is written according to the algorithm on
    P. 456 of Numerical Analysis 9th ED by Richard Burden
"""


from gaussian_elimination_2 import *

def norm(vec):
    return sqrt(dot(vec, vec))
    

def GS(A, b, x0, tolerance, max_iter):
    """
    not finished
    """


    n = rows(A)
    
    
    # Step 1
    k = 0
    error = [-1 for _ in range(max_iter)]    

    # Step 2
    while k < max_iter:
        x = [0 for _ in range(n)]

        # Step 3
        for i in range(n):
            x[i] = 1.0/A[i][i] * ( -sum( [A[i][j] * x[j] for j in range(i)])
                                 -sum( [A[i][j] * x0[j] for j in range(i+1, n)])
                                 +b[i])
            print("\nx")
            print x
            
        # Step 4
        error[k] = norm(subVectors(x, x0))
        if error[k] < tolerance:
            print("The procedure was successful")
            return [x, error]

        # Step 5
        k = k + 1

        # Step 6
        x0 = x

    print("Maximum number of iterations exceeded")
    return [x, error]

A = [[1,2,3],
     [3,4,10],
     [2,3,1]]

b = [1,2,3]
x0 = [0,0,0]
tol = 0.1
maxIter = 5
x = GS(A, b, x0, tol, maxIter)
