#! /usr/bin/python3.1
# -*- coding: utf-8 -*-

from math import sqrt,sin,cos,tan,pi,log,exp
from copy import copy, deepcopy
from functools import reduce
eps = 0.0001

def RungeKutta( f, g, segment, h, y0, z0):
	"""Iterative methods for the approximation of solutions of ODE"""
	x0 = segment[0]
	func = [ (x0,y0,z0) ]
	
	while(x0 < segment[1] - eps):
		k1 = h * f((x0, y0, z0))
		l1 = h * g((x0, y0, z0))
		k2 = h * f((x0 + h/2, y0 + k1/2, z0 + l1/2));
		l2 = h * g((x0 + h/2, y0 + k1/2, z0 + l1/2));
		k3 = h * f((x0 + h/2, y0 + k2/2, z0 + l2/2));
		l3 = h * g((x0 + h/2, y0 + k2/2, z0 + l2/2));
		k4 = h * f((x0 + h, y0 + k3, z0 + l3));
		l4 = h * g((x0 + h, y0 + k3, z0 + l3));
		dy = (k1 + 2 * k2 + 2 * k3 + k4) / 6;
		dz = (l1 + 2 * l2 + 2 * l3 + l4) / 6;
		x0 += h
		y0 += dy
		z0 += dz
		func.append((x0,y0,z0))
	return func

		
class Tridiagonal_Matrix:
	def __init__(self, f = ""):
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
		


def f(x):
	return x*x*x + 3*x + 1

def y(p):
	return p[1]


def f1(p):
	(x,y,z) = p
	return z

def f2(p):
	(x,y,z) = p
	return exp(x) + sin(y)

def pr(X):
	for xx in X :
		print("%.2f\t"%xx, end = '')
	print()

def shoting(f, g, segment, h, y0, y1 ):
	eta = []
	eta.append(1.0)
	eta.append(0.8)
	Rk0 = Rk1 = RungeKutta(f, g, (0, 1) , 0.1, 1, eta[-2]) [-1]
	
	delta1 = 1
	count = 0
	while (delta1 > eps)and(count < 9000):
		count+=1
		Rk0 = Rk1
		Rk1 = RungeKutta(f, g, segment , h, y0, eta[-1]) [-1]
		
		delta1 = abs( y(Rk1) - y1)
		delta0 = abs( y(Rk0) - y1)
		
		new_eta = eta[-1] - (y(Rk1)-y1)*(eta[-1] - eta[-2])/(y(Rk1) - y(Rk0))
		eta.append(new_eta)
		
		#~ print( "y(Rk0), y(Rk1)"  )
		#~ print( y(Rk0), y(Rk1)  )
	
	return RungeKutta(f1, f2, segment , h, y0, eta[-1])
	
	
def finite_diff(X, p, q, f, h, n = 5):
	ya = 1.0
	
	#~ May be +h*h !!!
	Tridiag = Tridiagonal_Matrix("")
	a = [ 0.0 ]
	b = [ -2 + h*h * q(X[1]) ]
	c = [ 1 - p(X[1])*h/2 ]
	d = [ h*h * f(X[1]) - 1 - p(X[1])*h/2.0*ya ]
	
	for k in range(2, n):
		a.append( 1 + p(X[k])*h/2.0 )
		b.append( -2 + h*h * q(X[k]) )
		c.append( 1 - p(X[k])*h/2 )
		d.append( h*h * f(X[k]) )
	
	#~ for first border
	#~ a.append( 1 - p(X[-1])*h/2.0 )
	#~ b.append( -2 + h*h * q(X[-1]) )
	#~ c.append( 0.0 )
	#~ d.append( h*h * f(X[k]) - 1 - p(X[-1])*h/2 *yb )
	
	a.append( 5)
	b.append(-7)
	c.append( 0)
	d.append( 0)
	
	#~ pr(a)
	#~ pr(b)
	#~ pr(c)
	#~ pr(d)
	
	
	Tridiag.a = a
	Tridiag.b = b
	Tridiag.c = c
	Tridiag.d = d
	Tridiag.n = len(a)
	
	
	y = [ya] + Tridiag.solve()
	return y
	
def RRR(f1, f2, k, p):
	# k = h1/h2, p - interpolation level
	return (f1 - f2) / (k**p - 1)
	

def main():

	p = lambda x: -x
	q = lambda x: -1
	f = lambda x: 0
	h = 0.2
	
	
	X = [ h*i for i in range(6)]
	X2 = [ h*i/2 for i in range(11)]
	X4 = [ h*i/4 for i in range(21)]
	print( "X = ", X)
	
	Sh = shoting(f1, f2, [0, 1] , 0.1, 1.0, 2.0)
	print("--Shoot--")
	print(Sh[-1])
	
	
	y1 = finite_diff(X, p, q, f, h, 5)
	y2 = finite_diff(X2, p, q, f, h/2, 10)
	y4 = finite_diff(X4, p, q, f, h/4, 20)
	
	
	print("---finite_diff---")
	print("X =")
	pr(X)
	print("y1 =")
	pr(y1)
	print("y2 =")
	pr(y2)
	print("y4 =")
	pr(y4[::2])
	
	
	print("eps:", RRR(y4[-1], y2[-1], 1/2, 2) )
	
#   y'' + p(x)y' + q(x) = 0
#	y(a) = y0; y(b) = y1;
	
	return 0

main()
