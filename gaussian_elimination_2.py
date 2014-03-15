"""
gauss_elim_1.py

Includes code for functions that do basic vector and
matrix arithmetic.  Most of these functions support
the function ge_1(aug) which takes an n by n+1
augmented matrix and returns a row-reduced matrix and
an approximate solution vector of the corresponding linear
system.  It uses gaussian elimination with a naive pivot
strategy.  That is, at each stage in the row reduction it
chooses, as the pivot, the first nonzero entry that lies
on or below the diagonal in the current column.

revision 0.01 02/06/12  added code [ds]
revision 0.02 02/06/12  gardened out some code [bic]
revision 1.23 02/08/12  added aug_2 example, cleaned code [ds]
revision 1.24 02/09/12  cleaned code some more [ds]
revision 2.00 03/12/12  changed filename, commented out tests [ds]

matrix examples: [[28, 12, 20, 28], [4, 32, 28, 16]] , 2 by 4
                 [[28, 12, 20, 28]] , 1 by 4
                 [[[28], [12], [20], [28]] , 4 by 1


vector example [28, 12, 20, 28] 

"""

# one way:


from math import sqrt
from decimal import *


def show(mat):
    "Print out matrix"
    outString = ""
##    for row in mat:
##        print(row)
    for row in mat:
        outString = outString + str(row) + "\n"
    print outString

def rows(mat):
    "return number of rows"
    return(len(mat))

def cols(mat):
    "return number of cols"
    return(len(mat[0]))
 
def zero(m,n):
    "Create zero matrix"
    new_mat = [[0 for col in range(n)] for row in range(m)]
    return new_mat
 
def transpose(mat):
    "return transpose of mat"
    new_mat = zero(cols(mat),rows(mat))
    for row in range(rows(mat)):
        for col in range(cols(mat)):
            new_mat[col][row] = mat[row][col]
    return(new_mat)

def dot(A,B):
    "vector dot product"
    if len(A) != len(B):
        print("dot: list lengths do not match")
        return()
    dot=0
    for i in range(len(A)):
        dot = dot + A[i]*B[i]
    return(dot)

def getCol(mat, col):
    "return column col from matrix mat"
    return([r[col] for r in mat])

def getRow(mat, row):
    "return row row from matrix mat"
    return(mat[row])

def setCol(mat, k, vec):
    "set kth column of matrix to be equal to vec"
    "operates in place"
    "vec must be a list of numbers, although the function does not check for this"
    for i in range(rows(mat)):
        mat[i][k] = vec[i]

def setRow(mat, k, vec):
    "seths kth row of matrix to be equal to vec"
    "operates in place"
    "vec must be a list of numbers of correct length"
    mat[k] = vec

def setDiagonal(mat, value):
    "sets each value on major diagonal to value"
    for row in range(rows(mat)):
        mat[row][row] = value

def identity(n):
    "creates nxn identity matrix"
    return [[0 if col != row else 1 for col in range(n)] for row in range(n)]

def matMult(mat1,mat2):
    "multiply two matrices"
    if cols(mat1) != rows(mat2):
        print("multiply: mismatched matrices")
        return()
    prod = zero(rows(mat1),cols(mat2))
    for row in range(rows(mat1)):
        for col in range(cols(mat2)):
            prod[row][col] = dot(mat1[row],getCol(mat2,col))
    return(prod)

def randIntMat(n, cust_seed = 1234, min = 0, max = 10):
    from random import randint
    from random import seed
    seed(cust_seed)
    return [[randint(min, max) for j in range(n)] for i in range(n)]
    

def vectorQ(V):
    "mild test to see if V is a vector"
    if type(V) != type([1]):
        return(False)
    if type(V[0]) == type([1]):
        return(False)
    return(True)

def scalarMult(a,mat):
    "multiply a scalar times a matrix"
    if vectorQ(mat):
        return([a*m for m in mat])
    for row in range(rows(mat)):
        for col in range(cols(mat)):
            mat[row][col] = a*mat[row][col]
    return(mat)

def addVectors(A,B):
    "add two vectors"
    if len(A) != len(B):
        print("addVectors: different lengths")
        return()
    return([A[i]+B[i] for i in range(len(A))])

def subVectors(A, B):
    "subtract two vectors"
    "Peter Varshavsky"
    if len(A) != len(B):
        print("subVectors: different lengths")
        return()
    return([A[i] - B[i] for i in range(len(A))])

### PV: why create a new matrix N?
def swaprows(M,i,j, s = None):
    "swap rows i and j in matrix M"
    "PV: added scale factor swapping" 
    N=copyMatrix(M)
    T = N[i]
    N[i] = N[j]
    N[j] = T

    # Swap scale factors
    if s != None:
        s[i], s[j] = s[j], s[i]
    
    return N

def copyMatrix(M):
    return([[M[row][col] for col in range(cols(M))]for row in
            range(rows(M))])

### PV: why create a new matrix N?
def addrows(M, f, t, scale=1):
    "add scale times row f to row t"
    N=copyMatrix(M)
    T=addVectors(scalarMult(scale,N[f]),N[t])
    N[t]=T
    return(N)

def mat2decimal(mat):
    "Convert matrix to decimal"
    return [[Decimal(col) for col in row] for row in mat]

### vectors vs rowVectors and colVectors
### the latter are matrices

def vec2rowVec(vec):
    "[a,b,c] -> [[a,b,c]]"
    return([vec])

def vec2colVec(vec):
    return(transpose(vec2rowVec(vec)))

def colVec2vec(mat):
    rowVec = transpose(mat)
    return(rowVec[0])

def augment(mat,vec):
    "given nxn mat and n length vector return augmented matrix"
    amat = []
    for row in range(rows(mat)):
        amat.append(mat[row]+[vec[row]])
    return(amat)


### The next two functions support checking a solution.

def getAandb(aug):
    "Returns the coef. matrix A and the vector b of Ax=b"
    m = rows(aug)
    n = cols(aug)
    A = zero(m,n-1)
    b = zero(m,1)
    for i in range(m):
        for j in range(n-1):
            A[i][j] = aug[i][j]
            
    for i in range(m):
        b[i] = aug[i][n-1]
    Aandb = [A,b]
    return(Aandb)

def checkSol_1(aug,x):
    "For aug=[A|b], returns Ax, b, and b-Ax as vectors"
    A  = getAandb(aug)[0]
    b  = getAandb(aug)[1]
    x_col_vec = vec2colVec(x)
    Ax = matMult(A,x_col_vec)
    r  = addVectors(b,scalarMult(-1.0,colVec2vec(Ax)))
    L  = [Ax,b,r]
    return(L)


### The naive gaussian elimination code begins here.

def findPivotrow1(mat,col):
    "Finds index of the first row with nonzero entry on or"
    "below diagonal.  If there isn't one return(-1)."
    for row in range(col, rows(mat)):
        if mat[row][col] != 0:
            return(row)
    return(-1)

def findPivotPP(mat, col):
    "Partial pivot by Peter Varshavsky"
    values = [mat[row][col] for row in range(col, rows(mat))]
    maxIndex = max(range(len(values)), key = lambda i: abs(values[i]))
    if abs(values[maxIndex]) > 0:
        return(maxIndex + col)
    else:
        return(-1)

def findPivotSPP(mat, col, s):
    "Scaled partial pivot by Peter Varshavsky"
    "unfinished"
    
    values = [Decimal(abs(mat[row][col]))/s[row] for row in range(col, rows(mat))]
    #print("Scale factors in SPP: ", s)
    #print("First column values in SPP: ", getCol(mat, 0))
    #print("Values in SPP: ", values)
    
    if max(values) > 0:
        return col + max(range(len(values)), key = lambda i: values[i])
    else:
        return(-1)
    
def scaleFactors(mat):
    "Scale factors by Peter Varshavsky"
    "Returns a vector of scale factors for matrix M"
    
    ncols = cols(mat)
    if rows(mat) == ncols:
        s = [abs(max(row, key = lambda x: abs(x))) for row in mat] # must change for absolute values
    else:
        s = [abs(max(row[0:-1], key = lambda x: abs(x))) for row in mat]
    return(s)


def rowReduce(M, pivotStrat = "naive"):
    "return row reduced version of M"
    "PV: added pivot strategy parameter"
    N = copyMatrix(M)
    cs = cols(M)-2   # no need to consider last two cols
    rs = rows(M)
    s = scaleFactors(M)
    #print("Scale factors: %s" %s)
    
    for col in range(cs+1):
        #print("Col: %s" %col)
        # find pivot
        if pivotStrat == "naive":
            #print "\nExecuting naive"
            j = findPivotrow1(N,col) # Naive Gaussian Elimination
        elif pivotStrat == "partial":
            #print "\nExecuting partial"
            j = findPivotPP(N, col) # Partial Pivoting
        elif pivotStrat == "scaled partial":
            #print "\nExecuting scaled partial"
            j = findPivotSPP(N, col, s)
        else:
            print("\nrowReduce: pivotStrat parameter value not accepted")

        # swap rows
        if j < 0:
            print("\nrowReduce: No pivot found for column index %d "%(col))
            return(N)
        else:
            if j != col:
                N = swaprows(N,col,j, s)
            scale = Decimal(-1.0) / N[col][col] # decimal
            for row in range(col+1, rs):
                N = addrows(N, col, row, scale * N[row][col])
                
    return(N)

def LUfactorize(M):
    "Author: Peter Varshavsky"
    "Returns a list [L, U]"
    "Has no checks for whether factorization is possible"

    cs = cols(M)-2
    rs = rows(M)

    L = identity(rs)
    U = copyMatrix(M)

    for col in range(cols(M)):
        atomic = identity(rs)
        atomic1 = identity(rs)
        
        for row in range(col+1, rs):
            "it seems redundant to use two atomic matrices"
            "one for knocking out columns in U"
            "and one to keep in L"
            atomic[row][col] = -U[row][col] / U[col][col]
            atomic1[row][col] = U[row][col] / U[col][col]
            
        U = matMult(atomic, U)
        L = matMult(L, atomic1)
    
    return([L, U])



    
def example2():
    print "************ EXAMPLE 2 FROM PAGE 404"
    A = [[1, 1, 0, 3],
         [2, 1,-1, 1],
         [3,-1,-1, 2],
         [-1,2, 3,-1]]
    L1 = [[1,0,0,0],
          [2,1,0,0],
          [3,4,1,0],
          [-1,-3,0,1]]

    
    U1 = [[1,1,0,3],
          [0,-1,-1,-5],
          [0,0,3,13],
          [0,0,0,-13]]
    
    #show(matMult(L1, U1))
    show(A)
    LU = LUfactorize(A)
    show(LU[0])
    show(LU[1])
    show(matMult(LU[0], LU[1]))
    print "*********** \EXAMPLE 2 FROM PAGE 404"

#example2()

"""
A = [[1, 1, 0, 3],
     [2, 1,-1, 1],
     [3,-1,-1, 2],
     [-1,2, 3,-1]]

L1 = [[1, 0, 0, 0],
      [-2, 1, 0, 0],
      [-3, 0, 1, 0],
      [1, 0, 0, 1]]

L2 = [[1, 0, 0, 0],
      [0, 1, 0, 0],
      [0, 1, 1, 0],
      [0, -2, 0, 1]]

L3 = [[1,0,0,0],
      [0,1,0,0],
      [0,0,1,0],
      [0,0,3,1]]
show(matMult(L3, A))
show(matMult(L1, L2))
print "grr"
show(matMult(matMult(L1, L2), L3))
show(matMult(L2, A))
show(matMult(L2, matMult(L1, A)))
"""

def backSub(M):
    """
    given a row reduced augmented matrix with nonzero 
    diagonal entries, returns a solution vector
    
    """
    cs = cols(M)-1 # cols not counting augmented col
    sol = [0 for i in range(cs)] # place for solution vector
    for i in range(1,cs+1):
        row = cs-i # work backwards
        sol[row] = ((M[row][cs] - sum([M[row][j]*sol[j] for
                    j in range(row+1,cs)])) / M[row][row]) 
    return(sol)


def diag_test(mat):
    """
    Returns True if no diagonal element is zero, False
    otherwise.
    
    """
    for row in range(rows(mat)):
        if mat[row][row] == 0:
            return(False)
    else:
        return(True)


def ge_1(aug, pivotStrat = "naive"):    
    """
    Given an augmented matrix it returns a list.  The [0]
    element is the row-reduced augmented matrix, and 
    ge_1(aug)[1] is a solution vector.  The solution
    vector is empty if there is no unique solution.
    
    """
    aug_n = rowReduce(aug, pivotStrat)
    if diag_test(aug_n):
        sol = backSub(aug_n)
    else:
        print("\nge_1(): There is no unique solution")
        sol = []
    results = [aug_n, sol]
    return(results)

def residueTest(mat, b, sol):
    "Error testing by Peter Varshavsky"
    "Returns a list with two components:"
    "a list of residues"

    prod = colVec2vec(matMult(mat, vec2colVec(sol))) # mat * sol
    
    residues = subVectors(b, prod)
    error = sqrt(dot(residues, residues))
    print "Solution: ", sol
    print "Error: %s" %error

    return error


### Some Testing code begins here.

A = [[4,-2,1],
    [-2,4,-2],
    [1,-2,4]]


C=[2,3,5]


aug = [[ 1.0, -1.0,  2.0, -1.0,  8.0],
       [ 0.0,  0.0, -1.0, -1.0, -4.0],
       [ 0.0,  2.0, -1.0,  1.0,  6.0],
       [ 0.0,  0.0,  2.0,  4.0, 12.0]]


aug_2 = [[ 1, -1,  2, -1,  8],
         [ 0,  0, -1, -1, -4],
         [ 0,  2, -1,  1,  6],
         [ 0,  0,  2,  2, 12]]



def showProcess(A,S):
    pivotStr = "scaled partial"
    print("\n***********************************************")
    print pivotStr
    print("***********************************************")
    
    "given matrix A and vector S, get B=AS and show solve for S"

    print("\n***********************************************")
    print("A")
    show(A)
    print("S=%s"%(S))
    B = matMult(A,vec2colVec(S))
    print("A*S=%s"%(B))
    AS = augment(A,colVec2vec(B))
    print("augment(A,A*S)=")
    show(AS)
    Ar=rowReduce(AS, pivotStr)
    print("row reduced Ar=")
    show(Ar)
    sol = backSub(Ar)
    print("solution is %s"% sol)
    
    print("\n***********************************************")
    print("aug = ")
    show(aug)
    aug_n = rowReduce(aug)
    print("aug_n = ")
    show(aug_n)
    print("the solution from ge_1(aug) is %s"%(ge_1(aug, pivotStr)[1]))
    print("\naug = ")
    show(aug)
    L = getAandb(aug)
    print("\nA = ")
    show(L[0])
    print("\nb = %s" %(L[1]))
    y = ge_1(aug, pivotStr)[1]
    print("\ny = ge_1(aug)[1] = %s"%(y))
    x_vec =  vec2rowVec(y)
    print("\nx_row_vec = %s"%(x_vec))
    x = transpose(vec2rowVec(y))
    print("\nx_col_vec = %s"%(x))
    Ax = matMult(L[0],x)
    print("\nAx = %s"%(Ax))
    print("\nAx_vec = %s"%(colVec2vec(Ax)))
    x_col_vec = vec2colVec(ge_1(aug, pivotStr)[1])
    L = checkSol_1(aug,ge_1(aug, pivotStr)[1])
    print("\n Ax = %s, \nb = %s, \nr = %s"%(L[0],L[1],L[2]))
    
    print("\n***********************************************")
    print("\naug_2 = ")
    show(aug_2)
    outputlist = ge_1(aug_2, pivotStr)
    print("\naug_2_n = ge_1(aug_2)[0] = \n")
    show(outputlist[0])
    print("\nge_1(aug_2)[1] is %s"%(outputlist[1]))

#showProcess(A,C)

### Peter testing things

    
def testDecimal():
    A = [[1, 1, 0, 3],
         [2, 1,-1, 1],
         [3,-1,-1, 2],
         [-1,2, 3,-1]]


    b2 = [1,2,3,4]
    Aaug = augment(A, b2)
    Adec = mat2decimal(Aaug)

    ans = ge_1(Adec)
    print "\nAnswer: \n",ans[1]


##    ans2 = ge_1(augA2, pivotStrat = "naive")
##    residueTest(A2, b2, ans2[1])
##    ans2a = ge_1(augA2, pivotStrat = "partial")
##    residueTest(A2, b2, ans2a[1])
##    ans2b = ge_1(augA2, pivotStrat = "scaled partial")
##    residueTest(A2, b2, ans2b[1])
##    print("Answer with naive", ans2[1])
##    print("Answer with PP", ans2a[1])
##    print("Answer with SPP", ans2b[1])
##    ans2wa = [-0.003125, 2.002, 2.002, 1.0]
##    residueTest(A2, b2, ans2wa)
    


def exercise2():
    A2 = [[10, 10, 10, 10**17],
          [1, 10**(-3), 10**(-3), 10**(-3)],
          [1, 1, 10**(-3), 10**(-3)],
          [1, 1, 1, 10**(-3)]]

    b2 = [10**17, 1, 2, 3]
    augA2 = augment(A2, b2)

    ans2 = ge_1(augA2, pivotStrat = "naive")
    residueTest(A2, b2, ans2[1])
    ans2a = ge_1(augA2, pivotStrat = "partial")
    residueTest(A2, b2, ans2a[1])
    ans2b = ge_1(augA2, pivotStrat = "scaled partial")
    residueTest(A2, b2, ans2b[1])
    print("Answer with naive", ans2[1])
    print("Answer with PP", ans2a[1])
    print("Answer with SPP", ans2b[1])
    ans2wa = [-0.003125, 2.002, 2.002, 1.0]
    residueTest(A2, b2, ans2wa)

#exercise2()

def decVec2float(vec):
    return [float(v) for v in vec]

def exercise2dec():
    A2 = [[10, 10, 10, 10**17],
          [1, 10**(-3), 10**(-3), 10**(-3)],
          [1, 1, 10**(-3), 10**(-3)],
          [1, 1, 1, 10**(-3)]]

    b2 = [10**17, 1, 2, 3]
    augA2 = augment(A2, b2)
    augA2dec = mat2decimal(augA2)
    
    ans2 = ge_1(augA2dec, pivotStrat = "scaled partial")
    print(decVec2float(ans2[1]))
    residueTest(A2, b2, decVec2float(ans2[1]))
##    residueTest(A2, b2, decVec2float(ans2[1]))
##    
##    ans2a = ge_1(augA2dec, pivotStrat = "partial")
##    residueTest(A2, b2, decVec2float(ans2a[1]))
##    
##    ans2b = ge_1(augA2dec, pivotStrat = "scaled partial")
##    residueTest(A2, b2, decVec2float(ans2b[1]))
##    
##    print("Answer with naive", ans2[1])
##    print("Answer with PP", ans2a[1])
##    print("Answer with SPP", ans2b[1])
##    #ans2wa = [-0.003125, 2.002, 2.002, 1.0]
    #residueTest(A2, b2, ans2wa)

getcontext().prec = 20
#exercise2dec()

def exercise2sympy():
    import sympy
    import numpy as np
    A2 = [[10, 10, 10, 10**17],
          [1, 10**(-3), 10**(-3), 10**(-3)],
          [1, 1, 10**(-3), 10**(-3)],
          [1, 1, 1, 10**(-3)]]

    b2 = [10**17, 1, 2, 3]
    augA2 = augment(A2, b2)
    
    A2sym = sympy.Matrix(augA2).rref()
    ansA2sym = [1.99609375, 0, 2, 1]
    residueTest(A2, b2, ansA2sym)

