#! /usr/bin/python3.1
# -*- coding: utf-8 -*-

from math import sqrt,sin,cos,atan,pi,log
from copy import copy, deepcopy
from matrix_2 import *

def f1(x):
	return (x[0][0]**2 + 4)*x[1][0] - 8

def f2(x):
	return (x[0][0]-1)**2 + (x[1][0]-1)**2 - 4

def phi1(x):
	return 8.0/(x[0][0]**2+4)

def phi2(x):
	return sqrt(4-(x[1][0]-1)**2)+1

def df1dx1(x):
	return 2*x[0][0]*x[1][0]
	
def df1dx2(x):
	return x[0][0]**2+4

def df2dx1(x):
	return 2*(x[0][0]-1)

def df2dx2(x):
	return 2*(x[1][0]-1)
	
	
def maximal(x):
	return max( x.transponate()[0] )
	
Af = [[df1dx1,df1dx2],[df2dx1,df2dx2]]
f = [f1,f2]

phis = [phi2,phi1]

segment = [[2.0,3.0],[0.5,1.0]]

m = 2
n = 2

def main():
	eps = 0.00001 #	eps = float(input("задайте eps: ") )
	dx = Matrix("",m,1)
	A = Matrix("",m,n)
	b = Matrix("",m,1)
	xk = Matrix("",m,1)
	xk.M = [ [(segment[i][0] + segment[i][1])/2.0] for i in range(m)]
		
	count = 0
	while 42:
		count+=1
		
		A.M = [ [ Af[i][j](xk)  for j in range(n)]
								for i in range(m)]
		
		b.M =[ [-f[i](xk)] for i in range(m)]
		
		A.b = b.transponate()[0]
		
		A.build_LU()
		v = A.solve(A.shift_b(A.b))
		
		dx.M = [ [v[i]] for i in range(dx.m)]
		xk = xk + dx
		
		if abs(maximal(dx))<eps:
			print ("ответ Ньютон:")
			print ("за %d шага" % count)
			xk.pr()
			break
			
	q = 0.5
	koef = q/(1-q)
	xk = Matrix("",m,1)
	xkm = Matrix("",m,1)
	xkm.M = [ [(segment[i][0] + segment[i][1])/2.0] for i in range(m)]

	count = 0
	while 42:
		count+=1
		xk.M = [ [phis[i](xkm)] for i in range(m)]
		if koef*abs((xk-xkm).norm()) < eps:
			print ("ответ Итерации:")
			print ("за %d шагов" % count)
			xk.pr()
			break
		xkm = deepcopy(xk)
			
	return 0

if (__name__ == "__main__"):
	main()
