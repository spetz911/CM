#!/usr/bin/env python3
# -*- coding: utf-8 -*-
###-------------------------------------------------------------------
### File    : elliptic.py
### Author  : Oleg Baskakov
### Description : elliptic pde solver
###
### 2011. Written for Moscow Aviation Institute.
###-------------------------------------------------------------------

from math import sqrt,sin,cos,tan,pi,log,exp,sqrt
from copy import copy, deepcopy
from functools import reduce
from tridiagonal import *

from pprint import pprint

from pde import *


from pprint import pprint
	
class Elliptic_PDE(PDE):

	approximate_boundary = '1lvl'
	method = 'liebmann'
	eps = 0.1
	tau = 0.0
	
	u_x = 0.0
	u = 0.0
	coef_x = [+1, -2, +1]
	coef_y = [+1, -2, +1]
	
	def __init__(self, pde = None):
		self.w = None
		self.max = 29
		MetaClass.copy(self, pde)
		MetaClass.print(self)
		self.grid = [[self.fun(x,y) for x in frange(0, self.lx, self.h)]
		                      for y in frange(0, self.ly, self.h)]

	def first_eq_2lvl(self, t):
		"""Find coefficients of first equation"""
		alpha = self.left[0]
		beta  = self.left[1]
		phi0  = self.left[2]
		U = self.grid[-1]
		
		a = self.u_xx
		b = self.u_x
		c = self.u
		h = self.h
		tau = self.tau
		
		a0 = 0
		b0 = alpha * (2*a*a/h + h/tau - c*h) - beta * (2*a*a - b*h)
		c0 = alpha * (-2*a*a/h)
		d0 = alpha * (U[0] * h/tau) - phi0(t) * (2*a*a - b*h)
		return (a0, b0, c0, d0)
	
	def last_eq_2lvl(self, t):
		"""Find coefficients of last equation"""
		alpha = self.right[0]
		beta  = self.right[1]
		phi1  = self.right[2]
		U = self.grid[-1]
		N = len(U)

		a = self.u_xx
		b = self.u_x
		c = self.u
		h = self.h
		tau = self.tau

		an = alpha * (-2*a*a/h)
		bn = alpha * (2*a*a/h + h/tau - c*h) + beta * (2*a*a + b*h)
		cn = 0
		dn = alpha * (U[-1] * h/tau) + phi1(t) * (2*a*a + b*h)
		return (an, bn, cn, dn)

	def boundary_eq(self, *args):
		if self.approximate_boundary == '1lvl':
			return self.boundary_eq_1lvl(*args)
		if self.approximate_boundary == '2lvl':
			return self.boundary_eq_2lvl(*args)

	def boundary_eq_1lvl(self, u, b_cond, tmp0, tmp1, tmp2):
		"""For x/y independency"""
		alpha = -b_cond[0]
		beta  = b_cond[1]
		phi   = b_cond[2]
		N = len(u)
		h = self.h

		koef = (beta - alpha/h)
		res = [0.0] + [(phi(i*h) - u[i]*alpha/h) / koef for i in range(1, N-1)] + [0.0]
	#	print("boundary")
	#	print_vec(u)
	#	print_vec([0.0] + [phi(i*h) for i in range(1, N-1)] + [0.0])
	#	print_vec(res)
		return res
	
	def boundary_eq_2lvl(self, u, b_cond, a, b, f):
		"""For x/y independency"""
		alpha = b_cond[0]
		beta  = b_cond[1]
		phi   = b_cond[2]
		N = len(u)
		h = self.h
		coef_x = self.coef_x
		
		res = [0.0]*N
		res[1:-1] = self.threads_pool.map(
			 lambda i, u: (phi(i*h) + alpha*u[i] +
		                   alpha/(a*h) * (b*scalar(coef_x, u[i-1:i+2]) - alpha*f(i*h))) /
		                  (alpha + beta), [(i,u) for i in range(1, N-1)])
		return res







	def Liebmann_method(self, method = "Liebmann"):
		"""Just solve equation"""
		U = deepcopy(self.grid[-1])
		N = len(U)
		M = len(U[0])
		h = self.h
		lx = self.lx
		ly = self.ly
		fun = self.fun
		
		##correct size
		for i in range(1, N-1):
			for j in range(1, M-1):
				U[i][j] = 0.25*(U[i][j-1] + U[i-1][j] + U[i+1][j] + U[i][j+1] - h*h * fun(i*h, j*h))
		
		#solve boundary conditions here!
#		print_mat(U)
		
		U[0] = self.boundary_eq(U[1], self.south, self.u_yy, self.u_xx, lambda z:fun(z,0))
		U[-1] = self.boundary_eq(U[-2], self.north, self.u_yy, self.u_xx, lambda z:fun(z,ly))
		U = [list(x) for x in zip(*U)] #transponate
		U[0] = self.boundary_eq(U[1], self.west, self.u_xx, self.u_yy, lambda z:fun(0,z))
		U[-1] = self.boundary_eq(U[-2], self.east, self.u_xx, self.u_yy, lambda z:fun(lx,z))
		U = [list(x) for x in zip(*U)] #transponate

		return U


	def Zeidel_method(self):
		"""Just solve equation"""
		U = deepcopy(self.grid[-1])
		N = self.N
		M = self.M
		h = self.h
		lx = self.lx
		ly = self.ly
		fun = self.fun

		tmp = [U[-2][j] for j in range(M)] 
		top = self.boundary_eq(tmp, self.north, self.u_yy, self.u_xx, lambda z:fun(z,ly))
		tmp = [U[i][1] for i in range(N)]
		left = self.boundary_eq(tmp, self.west, self.u_xx, self.u_yy, lambda z:fun(0,z))
		for j in range(1, M-1):
			U[-1][j] = top[j]
		for i in range(1, N-1):
			U[i][0] = left[i]
		
		# U[i][j-1], U[i-1][j] takes from this iteration
		# U[i+1][j], U[i][j+1] takes from last iteration
		for j in range(1, M-1):
			for i in range(1, N-1):
				U[i][j] = 0.25*(U[i][j-1] + U[i-1][j] + U[i+1][j] + U[i][j+1] - h*h * fun(i*h, j*h))
		
		tmp = [U[1][j] for j in range(M)] 
		footer = self.boundary_eq(tmp, self.south, self.u_yy, self.u_xx, lambda z:fun(z,0))
		tmp = [U[i][-2] for i in range(N)]
		right = self.boundary_eq(tmp, self.east, self.u_xx, self.u_yy, lambda z:fun(0,z))
		for j in range(1, M-1):
			U[0][j] = footer[j]
		for i in range(1, N-1):
			U[i][-1] = right[i]

		return U

	def SOR(self):
		"""Successive over-relaxation"""
#		print('SOR')
		w = self.w
		U0 = self.grid[-1]
		U1 = self.grid[-2]
		n = len(U0)
		m = len(U0[0])
		
		for i in range(n):
			for j in range(m):
				U0[i][j] = U1[i][j] + w * (U0[i][j] - U1[i][j])


	def estimate_error(self):
		U0 = self.grid[-1]
		U1 = self.grid[-2]
		
		eps = max(max(abs(e0-e1) for e0,e1 in zip(v0,v1))
		                         for v0,v1 in zip(U0,U1))
		return eps
	
	def solve(self, method = 'liebmann'):
		"""Description"""
#		self.M = int((self.lx+0.00001) / self.h)
#		self.N = int((self.ly+0.00001) / self.h)

		self.M = len(list(frange(0, self.lx+0.00001, self.h)))
		self.N = len(list(frange(0, self.ly+0.00001, self.h)))

		self.grid = [[0.0]*self.M for i in range(self.N)]
		Us = [self.grid]
		self.grid = Us
		self.count = 1
		
		if method == 'liebmann':
			one_iteration = self.Liebmann_method
		elif method == 'zeidel':
			one_iteration = self.Zeidel_method
		elif method == 'relax':
			if not self.w:
				self.w = 1.5
			one_iteration = self.Liebmann_method
		else:
			print("method %s not supported" % method)
			return
		
		Us.append(one_iteration())
		while self.estimate_error() > self.eps and self.count < self.max:
			
			Us.append(one_iteration())
			if self.w: self.SOR()
#			print_mat(Us[-1])
			self.count += 1

		print("complete in %d iters" % (len(Us)-1))
		print("estimate", self.estimate_error())
		return Us



#=====================================================================

def main():
	pass

if __name__ == "__main__":
	main()
