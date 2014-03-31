import numpy as np
import matplotlib.pyplot as plt
from math import *

X = [0.00, 1/6.0, 1/3.0, 1/2.0, 7/12.0, 2/3.0, 3/4.0, 5/6.0, 11/12.0, 1.00]
dx = 1/6.0

def Lagrange_Interpolation(X, dx):
    x1 = [0.0]
    y = [(1.6 * pow(e, -2 * i)) * (sin(3 * pi * i)) for i in X]
    
    for i in range(6):
        element = x1[i]
        add = element + dx
        x1.append(add)

    y1 = [(1.6 * pow(e, -2 * i)) * (sin(3 * pi * i)) for i in x1]

    plt.plot(X, y, 'ro-', x1, y1, 'bs-')
    plt.axis([0, 1, -1, 2])
    plt.xlabel("X values")
    plt.ylabel("Y values")
    plt.title('Lagrange Polynomials')
    plt.text(0.26, 1.7, r'$y = 1.6exp(-2x)sin(3 \pi x)$', fontsize = 20)
    plt.text(0.25, 1.25, r'$Blue Graph: Lagrange Polynomial$', fontsize = 15)
    plt.grid(True)         
        
    plt.show()

Lagrange_Interpolation(X, dx)
