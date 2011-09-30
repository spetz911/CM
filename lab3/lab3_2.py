#! /usr/bin/python3.1
# -*- coding: utf-8 -*-

#~ import matplotlib.pyplot as plt
#~ import numpy as np

from math import sqrt,sin,cos,atan,pi,log
from copy import copy, deepcopy
from functools import reduce

class Tridiagonal_Matrix:
	def __init__(self):
			self.a = []
			self.b = []
			self.c = []
			self.d = []
			self.n = 0

	def solve(self):
		"""Method progonki"""
		
		a = self.a
		b = self.b
		c = self.c
		d = self.d
		n = len(d)
		
		P = []
		Q = []
		P.append(-c[0]/b[0])
		Q.append( d[0]/b[0])
		
		for i in range(1, n):
			P.append( -c[i] / (b[i]+a[i]*P[i-1]) ) 
			Q.append( (d[i] - a[i]*Q[i-1]) / (b[i] + a[i]*P[i-1]) )

		x = [0]*n
		x[n-1] = Q[n-1]
		for i in range(n-2, -1, -1):
			x[i] = P[i]*x[i+1] + Q[i]
			
		return x

Xi = [0.0, 1.0, 2.0, 3.0, 4.0]
#Yi = [0.0, 0.5, 0.86603, 1.0, 0.86603]
Yi = [0.0, 1.8415, 2.9093, 3.1411, 3.2432]

Hi = [0.0] + [(Xi[i] - Xi[i-1]) for i in range(1,len(Xi))]

X0 = 1.5
Y0 = 2.4969

def roundx(v):
	return [round(xx, 5) for xx in v]


def build_ci(Y, X, H):
	n = len(X)
	M = Tridiagonal_Matrix()	
	
	for i in range(2, n): # numerating from zero
		a = H[i-1]
		b = 2*(H[i-1]+H[i])
		c = H[i]
		d = 3*( (Y[i]-Y[i-1])/H[i] - (Y[i-1]-Y[i-2])/H[i-1] )
		if (i==2)  : a = 0
		if (i==n-1): c = 0
		M.a.append(a)
		M.b.append(b)
		M.c.append(c)
		M.d.append(d)
	
	print("Tridiagonal_Matrix:")
	print(roundx(M.a))
	print(roundx(M.b))
	print(roundx(M.c))
	print(roundx(M.d))
	
	x = [0.0] + M.solve()
	
	print ("x = ", roundx(x) )
	
	return x
	
	
def spline1(A, B, C, D, X, x):
	""" calculate spline"""
	if (x < X[0]) : return 0.0
	if (x > X[-1]) : return 0.0
	
	i = 0
	while (X[i+1] < x): i+=1
	# segment X[i]..X[i+1]
	x1 = X[i]
	return A[i] + B[i]*(x-x1) + C[i]*(x-x1)**2 + D[i]*(x-x1)**3
		
def cr_spline(A, B, C, D, X):
	""" create spline"""
	return lambda x: spline1(A, B, C, D, X, x)




def main():
	print("table function: ")
	print(Xi)
	print(Yi)

	Ai = Yi[0 : -1]
	Ci = build_ci(Yi, Xi, Hi)
	
	n = len(Ci)-1
	Bi = [( (Yi[i+1]-Yi[i])/Hi[i+1] - Hi[i+1]*(Ci[i+1]+2*Ci[i])/3 ) for i in range(0,n)]
	Bi.append( (Yi[n+1]-Yi[n])/Hi[n] - Hi[n]*Ci[n]*2/3 )
	
	Di = [( (Ci[i+1]-Ci[i])/(3*Hi[i+1]) ) for i in range(0,n)]
	Di.append( -Ci[n]/(3*Hi[n]) )
	
	print("proverka:")
	print("A = ", roundx(Ai))
	print("B = ", roundx(Bi))
	print("C = ", roundx(Ci))
	print("D = ", roundx(Di))
	
	sp = cr_spline(Ai, Bi, Ci, Di, Xi)
	
		
	#~ X = np.arange(Xi[0], Xi[-1], 0.01)
	#~ Y = [ sp(xx) for xx in X]

	#~ plt.plot(X,Y,color="red")
	#~ plt.plot(Xi,Yi,color="grey")
	#~ plt.show()
		
	
	print("\ntest1 ", Y0,"=", sp(X0) )
	
	return 0

main()
