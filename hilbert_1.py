"""
hilbert_1.py

Includes code for functions that support basic vector and
matrix arithmetic.  It also includes functions that generate
Hilbert Matrices, and the inhomogeneous b-vector needed in
Assignment 1.

    Revision 1.00 03/12/13 ds

"""
from math import *

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


def show(mat):
    "Print out matrix"
    for row in mat:
        print(row)

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


### Hilbert matrix and related functions

def hilbert_b(n):
    # creates "Hilbert" b vector
    b_vec = [0 for k in range(n)]
    b_vec[0] = 2/pi
    b_vec[1] = 1/pi
    for k in range(2,n):
        b_vec[k] = 1/pi -(k*(k-1)/pi**2)*b_vec[k-2]
    return(b_vec)

def hilbert(n):
    # creates a Hilbert matrix
    h_n = zero(n,n)
    for row in range(n):
        for col in range(n):
            h_n[col][row] = 1.0/(1.0+row+col)
    return(h_n)

def n_C_k(n,k):
    # computes the binomial coefficien nCk
    if k > n:
        n_ch_k = 0
    prod = 1.0
    for j in range(k):
        prod = prod*(1.0*(n-j)/(k-j))
    return(prod)

def hilbert_inv(n):
    # creates the inverse of a Hilbert matrix
    h_inv = zero(n,n)
    for k in range(n):
        for m in range(n):
            h_inv[k][m] = ((-1)**(k+m))*((k+m+1)*
                                          (n_C_k(n+k,n-m-1))*
                                          (n_C_k(n+m,n-k-1))*
                                          (n_C_k(k+m,k))**2 )
    return(h_inv)




