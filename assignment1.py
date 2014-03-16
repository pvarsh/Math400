"""
assignment_1_group_8.py

3/15/14  bic M400

"""

#######  START Administrivia 
m400group = 8   # change this to your group number

m400names = ['Deborah Castro', 'Yueming Liu', 'Mirai Furukawa', 'Peter Varshavsky'] # change this for your names

def printNames():
    print("assignment_1_group_8.py for group %s:"%(m400group)),
    for name in m400names:
        print("%s, "%(name)),
    print

printNames()

#######  END Administrivia



import matplotlib.pyplot as plt
from numpy import arange
from math import *

#from decimal import *


def show(mat):
    "Print out matrix"
    outString = ""
    
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
    "set kth row of matrix to be equal to vec"
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
    "Author: Peter Varshavsky"
    "generate an nxn matrix of pseudorandom integers from min to max"
    
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

def scalarMult(a, mat):
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
    "create a copy of a matrix"
    return([[M[row][col] for col in range(cols(M))]for row in
            range(rows(M))])

### PV: why create a new matrix N?
def addrows(M, f, t, scale=1):
    "add scale times row f to row t"
    N=copyMatrix(M)
    T=addVectors(scalarMult(scale,N[f]),N[t])
    N[t]=T
    return(N)
"""
def mat2decimal(mat):
    "Convert matrix to decimal"
    return [[Decimal(col) for col in row] for row in mat]
"""
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

def augmentMat(mat1, mat2):
    amat = []
    for i in range(rows(mat1)):
        amat.append(mat1[i] + mat2[i])
    return amat
        

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

def L2norm(vec):
    return sqrt(dot(vec, vec))

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

def findPivotrow1(mat, col):
    "Finds index of the first row with nonzero entry on or"
    "below diagonal.  If there isn't one return(-1)."
    "Used in Naive Gaussian Elimination"
    for row in range(col, rows(mat)):
        if mat[row][col] != 0:
            return(row)
    return(-1)

def findPivotPP(mat, col):
    "Partial pivot by Peter Varshavsky"
    "Used in Gaussian Elimination with Partial Pivoting"
    values = [mat[row][col] for row in range(col, rows(mat))]
    maxIndex = max(range(len(values)), key = lambda i: abs(values[i]))
    if abs(values[maxIndex]) > 0:
        return(maxIndex + col)
    else:
        return(-1)

def findPivotSPP(mat, col, s):
    "Scaled partial pivot by Peter Varshavsky"
    
    #values = [Decimal(abs(mat[row][col]))/s[row] for row in range(col, rows(mat))]
    values = [float(abs(mat[row][col]))/s[row] for row in range(col, rows(mat))]


    
    if max(values) > 0:
        return col + max(range(len(values)), key = lambda i: values[i])
    else:
        return(-1)
    
def scaleFactors(mat):
    "Scale factors by Peter Varshavsky"
    "Returns a vector of scale factors for matrix M"
    "Used in Scaled Partial Pivot"
    
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
                
            #scale = Decimal(-1.0) / N[col][col] # decimal
            scale = -1.0 / N[col][col] # decimal
            for row in range(col+1, rs):
                N = addrows(N, col, row, scale * N[row][col])
                
    return(N)

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


def forwardSub(M):
    """
    given a row reduced augmented lower-triangular matrix with nonzero
    diagonal entries, returns a solution vector
    """
    cs = cols(M) - 1
    sol = [0 for i in range(cs)]
    for row in range(0, cs):
        sol[row] =   (float(M[row][cs]) - sum( [M[row][j] * sol[j] for j in range(row)] )) / M[row][row]

    return(sol)

def solveLU(A, b):
    L, U = LUFact(A)
    y = forwardSub(augment(L, b))
    x = backSub(augment(U, y))

    return x



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

def gauss_seidel(A, b, x0, tolerance, max_iter):
    """
    Author: Peter Varshavsky
    """
    
    n = rows(A)
    
    # Step 1
    k = 0
    error = [-1 for _ in range(max_iter + 1)]    

    # Step 2
    while k <= max_iter:
        x = [0 for _ in range(n)]

        # Step 3
        for i in range(n):
            #print "i: ", range(i)
            #print "i+1 ... n: ", range(i+1, n)
            x[i] = 1.0/A[i][i] * ( -sum( [A[i][j] * x[j] for j in range(i)])
                                   -sum( [A[i][j] * x0[j] for j in range(i+1, n)])
                                   +b[i] )
        #print("\nx%s :" %k)
        #print x
        #print("\nError: ")
        #print error
            
        # Step 4
        error[k] = L2norm(subVectors(x, x0))
        if error[k] < tolerance:
            print("The procedure was successful")
            return [x, error]

        # Step 5
        k = k + 1

        # Step 6
        x0 = x

    print("Maximum number of iterations exceeded")
    return [x, error]

def LUFact(M):
    "Author: Lauren"
    "return row reduced version of L"
    N = copyMatrix(M)
    L = zero(len(M),len(M[0]))
    cs = cols(M)-2   # no need to consider last two cols
    rs = rows(M)

    #create identity matrix out of L
    for i in range(len(L)):
        for j in range(len(L)):
            if i == j:
                L[i][j] = 1

    for col in range(cs+1):
        scale = -1.0 / N[col][col]
        for row in range(col+1,rs):
            L[row][col] = (-1)*scale*N[row][col]
            N=addrows(N, col, row, scale * N[row][col])
    for col in range(cs+1):
        for row in range(col+1,rs):
            if N[row][col]<10**-14:
                N[row][col]=0         

    """
    print"L ="
    show(L)
    print "U ="
    show(N)
    print""
    """
    return(L,N)

def residueTest(mat, b, sol, verbose = True):
    "Error testing by Peter Varshavsky"
    "Returns a list with two components:"
    "a list of residues"

    prod = colVec2vec(matMult(mat, vec2colVec(sol))) # mat * sol
    
    residues = subVectors(b, prod)
    error = sqrt(dot(residues, residues))

    if verbose:
        print "Residues: ", residues
        print "Solution: ", sol
        print "Error: %s" %error

    return error


def makeHilbert(n):
    " Returns an nxn Hilbert matrix "
    " Author: PV "
    return [[1.0/(col+1) for col in range(row, n+row)] for row in range(0, n)]

def makeCvector(n):
    " Returns a vector of $\int_0^1 x^k \sin(\pi x)$ for k from 0 to n "
    " Author: PV "
    from math import pi
    
    b = [2/pi]
    b.append(1/pi)

    for i in range(2, n):
        b.append(1/pi - i*(i-1)*b[-2]/(pi**2))
    return b




def exercise2():
    print "\n****************************\nEXERCISE 2\n****************************"

    exponent = 17.45
    
    A = [[10, 10, 10, 10**exponent],
         [1, 10**(-3), 10**(-3), 10**(-3)],
         [1, 1, 10**(-3), 10**(-3)],
         [1, 1, 1, 10**(-3)]]

    b = [10**exponent, 1, 2, 3]
    augA = augment(A, b)

    print "This exercise compares the solutions of a system of equations Ax = b"
    print "using Gaussian Elimination with three pivoting strategies\n"
    
    print "Matrix A: "    
    show(A)
    print "Vector b: "
    print(b)

    print "\nNaive Gaussian Elimination"
    ans_naive = ge_1(augA, pivotStrat = "naive")
    residueTest(A, b, ans_naive[1])

    print "\nGaussian Elimination with Partial Pivoting"
    ans_partial = ge_1(augA, pivotStrat = "partial")
    residueTest(A, b, ans_partial[1])
    

    print "\nGaussian Elimination with Scaled Partial Pivoting"
    ans_spp = ge_1(augA, pivotStrat = "scaled partial")
    residueTest(A, b, ans_spp[1])

    
    
#exercise2()

def exercise3():
    print "\n****************************\nEXERCISE 3\n****************************"

    solGE = []
    solGESPP = []
    solLU = []
    solGS = []


    maxN = 25
    for n in range(1, maxN + 1):
        H = makeHilbert(n)
        b = makeCvector(n)

        #Gaussian Elimination
        _, sol = ge_1(augment(H, b), pivotStrat = 'naive')
        error = residueTest(H, b, sol, verbose = False)
        solGE.append([sol, error])
        #print "\nGE solution:"
        #print sol

        #Gaussian Elimination with Scaled Partial Pivoting
        _, sol = ge_1(augment(H, b), pivotStrat = "scaled partial")
        error = residueTest(H, b, sol, verbose = False)
        solGESPP.append([sol, error])
        #print "\nGESPP solution:"
        #print sol

        #Gauss-Seidel
        x0 = [0 for _ in range(n)]
        tolerance = 0.1
        max_iter = 10000
        sol = gauss_seidel(H, b, x0, tolerance, max_iter)[0]
        error = residueTest(H, b, sol, verbose = False)
        solGS.append([sol, error])
        #print "\nGS solution:"
        #print sol

        #LU factorization
        sol = solveLU(H, b)
        error = residueTest(H, b, sol, verbose = False)
        solLU.append([sol, error])
        #print "\nLU solution:"
        #print sol

   
    plt.plot( arange(1, maxN + 1, 1), getCol(solGE, 1), color = 'black', label = 'Naive GE')
    plt.plot( arange(1, maxN + 1, 1), getCol(solGS, 1), color = 'orange', label = 'Gauss-Seidel')
    plt.plot( arange(1, maxN + 1, 1), getCol(solGESPP, 1), color = 'red', label = 'SPP GE')
    plt.plot( arange(1, maxN + 1, 1), getCol(solLU, 1), color = 'cyan', label = 'LU')
    plt.legend(loc='upper left')
    plt.axis([0, maxN , 0, 0.0000004], 'equal')

    plt.show()

    #for i in [1,5,10,15,20,25]:
        
    """
    print "\nNaive Gaussian Elimination errors:"
    print getCol(solGE, 1)
    print "\nGaussian Elimination with Scaled Partial Pivoting errors:"
    print getCol(solGESPP, 1)
    print "\nGauss-Seidel errors:"
    print getCol(solGS, 1)
    """



#exercise3()

def exercise4():
    print "\n****************************\nEXERCISE 4\n****************************"

    n = 4
    H = makeHilbert(n)
    I = identity(n)

    HI = augmentMat(H, I)
    show(HI)

    rowReduce(HI)
    #show(I)
    
    #show(H)
    
    
exercise4()
