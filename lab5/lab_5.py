#! /usr/bin/python3.1
# -*- coding: utf-8 -*-

from math import sqrt,sin,cos,tan,pi,log,exp
from copy import copy, deepcopy
from functools import reduce
eps = 0.0001

def frange(x0, x1, d):
	return [ x0 + d*k for k in range(int((x1-x0+eps/10)/d)+1)]

def print_mat(mat):
	print("matrix:")
	for row in mat:
		for elem in row:
			print(round(elem, 3), end = "\t")
		print()
	print("--------")

def print_vec(row):
	print("vector:")
	for elem in row:
		print(round(elem, 4), end = "\t")
	print()
	print("--------")

class Tridiagonal_Matrix:
	def solve(self):
		"""Method progonki"""
		a = self.a
		b = self.b
		c = self.c
		d = self.d
		n = self.n
		
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

def explicit_method( u, step):
	N = len(u)
	sigma = tau * alpha**2 / h**2
	return [border1[step*tau]] +
	       [sigma*(u[i-1] - 2*u[i] + u[i+1]) for i in range(1, N-1)]+
	       [border2[step*tau]]

def implicit_method():
	sigma = alpha**2 * tau / h**2
	N = len(u)
	
	a = sigma
	b = -(1+2*sigma)
	c = sigma
	d = -u[j]
	
	M = Tridiagonal_Matrix();
	
	M.a = [a] * N
	M.b = [b] * N
	M.c = [c] * N
	M.d = [d] * N
	M.n = len(N)
	
	M.a[1] = 0
	M.c[N] = 0
	M.d[1] = -(u[1] + sigma*border1((step+1)*tau))
	M.d[N] = -(u[N] + sigma*border2((step+1)*tau))
	
	print_mat([M.a, M.b, M.c, M.d])

	x = M.solve()
	print_vec( [u0] + x + [u1])

start_f = lambda x: sin(x)
border1 = lambda x: 0
border2 = lambda x: 0

x0 = 0
x1 = 1

h = 0.1
tau = 0.5

T = 7.0

u0 = [start_f(x) for x in frange(x0, x1, h)]
U1 = [u0]
U2 = [u0]

for (step,t) in enumerate(frange(0, T, tau)):
	U1.append( implicid_method(U[-1], step))
	U2.append( explicid_method(U[-1], step))



#=====================================================================

