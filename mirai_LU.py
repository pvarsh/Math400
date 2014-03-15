

from gaussian_elimination_2 import *


def lUFactorization(mat):
    "return L[0] U[1]"
    N = copyMatrix(mat)
    cs = cols(mat)
    rs = rows(mat)
    L = []
    for col in range(cs):
        rowVec = [0.0 for i in range(col)]
        scale = -1.0 / N[col][col]
        for j in range(col,rs):
            rowVec.append(N[j][col]/N[col][col])
        for row in range(col+1,rs):
            if col != cs:
                N=addrows(N, col, row, scale * N[row][col])
        L.append(rowVec)
    return(transpose(L),N)

A = [[1,1,0,3],
     [2,1,-1,1],
     [3,-1,-1,2],
     [-1,2,3,-1]]
lu = lUFactorization(A)

show(lu[0])
show(lu[1])


# Making a change to show AP

print "testing randIntMat"
show(randIntMat(n = 5, cust_seed = 123, min = 0, max = 4))
