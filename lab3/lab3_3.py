#! /usr/bin/python3.1
# -*- coding: utf-8 -*-

from math import sqrt,sin,cos,atan,pi,log
from copy import copy, deepcopy
from functools import reduce
from matrix_2 import *

import matplotlib.pyplot as plt
import numpy as np


#Xi = [-1.0, 0.0, 1.0, 2.0    , 3.0, 4.0    ]
#Yi = [-0.5, 0.0, 0.5, 0.86603, 1.0, 0.86603]

Xi = [0.0, 1.7   , 3.4   , 5.1   , 6.8   , 8.5   ]
Yi = [0.0, 1.3038, 1.8439, 2.2583, 2.6077, 2.9155]

X0 = 1.5
Y0 = 2.4969

def roundx(v):
	return [round(xx, 4) for xx in v]

def P_n(A):
	return lambda x: sum([ A[i]*x**i for i in range(0,len(A)) ] )


def cr_MNK(Y, X, m = 2):
	n = m + 1

	A = Matrix("", n, n)
	A.b = [sum(Y)] + [ sum([(y_j * x_j**k) for (y_j, x_j) in zip(Y, X)])
												for k in range(1,n)    ]
	
	A.M = [ [sum([ x_j**(i+k) for x_j in X ])
								for k in range(n) ]
									for i in range(n) ]

#	A.M[0][0] = len(X) #~ does not matter
	A.pr()
	
	print("b = ", roundx(A.b) )
	
	A.n = len(A.M)
	A.m = len(A.M)
	A.build_LU()
	
	print("p = ", A.p)
	
	Ai = A.solve(A.shift_b(A.b))
	
	print("Ai = ", roundx(Ai) )
	return Ai
	
def F_eps(f, Y, X):
	return sum( [(f(xx)-yy)**2 for (xx,yy) in zip(X, Y)] )
 

def main():
	print("table function: ")
	print(Xi)
	print(Yi)

	
	Ai = cr_MNK(Yi, Xi, 2)
	polynom = P_n(Ai)


	X = np.arange(Xi[0], Xi[-1], 0.01)
	Y = [ polynom(xx) for xx in X]

	plt.plot(X,Y,color="red")
	plt.plot(Xi,Yi,color="grey")
	plt.show()

		



	print("\ntest1 ", Yi[1],"=", polynom(Xi[0]) )
	print("F_bolshoe     =", F_eps(polynom, Yi, Xi) )
	
	return 0

main()
