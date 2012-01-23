#! /usr/bin/python3.1
# -*- coding: utf-8 -*-

from math import sqrt,sin,cos,tan,pi,log,exp
from copy import copy, deepcopy
from functools import reduce
eps = 0.0001

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

def el_func(x, el, grid):
	if not( 0 <= el <= len(grid)): return 0
	
	z = -abs(x - grid[el])/(grid[1] - grid[0]) + 1.0
	if z<0: z = 0
	return z

class my_fun:
	def __init__(self, el = 0, grid = [0.0,1.0], f = None):
		if f:
			self.f = f
		else:
			self.f =  lambda x: el_func(x, el, grid)
	
	def D(self):
		eps = 0.0001
		fun = lambda x: (self.f(x+eps) - self.f(x-eps))/(2*eps)
		return my_fun(f = fun)
	
	def __neg__(self):
		fun = lambda x: -self.f(x)
		return my_fun(f = fun)

	def __add__(self, func2):
		fun = lambda x: self.f(x) + func2.f(x)
		return my_fun(f = fun)
		
	def __sub__(self, func2):
		fun = lambda x: self.f(x) - func2.f(x)
		return my_fun(f = fun)

	def __mul__(self, func2):
		fun = lambda x: self.f(x) * func2.f(x)
		return my_fun(f = fun)

	def integral(self, seg):
		d = eps
		func = self.f
		res = sum([func(x)*d for x in frange(seg[0], seg[1], d)])
		return res

seg = (0.0, 1.0)
n = 3

grid = [ seg[0] + k*(seg[1]-seg[0])/n for k in range(n+1) ] 


Func = [ my_fun(el, grid) for el in range(n+1)]


# print( integral(df(Func[0]), (0.0,0.7)) )

def matrix_elem(F, i, j):
	return ((my_fun.D(F[i]) * my_fun.D(F[j])) + (F[i] * F[j]))


def matrix_k( Func, grid, el):
	x0 = (grid[el] + grid[el])/2
	M = [[my_fun.integral(matrix_elem(Func,i,j), (grid[el],grid[el+1]))
					for j in [el, el+1] ]
					for i in [el, el+1] ]
	return M


system = [[ 0.0 for i in range(n+1)]
                for j in range(n+1)]

for k in range(n):
	mat = matrix_k(Func, grid, k)
	for i in [0,1]:
		for j in [0,1]:
			system[i+k][j+k] += mat[i][j]

print_mat(system)

M = Tridiagonal_Matrix();
a = [0.0]+[ system[i][i-1] for i in range(1,n+1)]
b = [ system[i][i]   for i in range(n+1)]
c = [ system[i][i+1] for i in range(0,n)]+[0.0]

print_mat([a,b,c])

u0 = 0.0
u1 = 1.0

# trunc first and last equation
M.a = a[1:-1]
M.b = b[1:-1]
M.c = c[1:-1]
M.d = [0.0] * len(M.a)
M.d[0] = -system[1][0] * u0
M.d[-1] = -system[-1][2] * u1
M.n = len(M.a)

print_mat([M.a,M.b,M.c, M.d])

x = M.solve()
print_vec( [u0] + x + [u1])

#=====================================================================

