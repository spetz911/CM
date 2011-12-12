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
	
	u_x = 0.0
	u = 0.0
	coef_x = [+1, -2, +1]
	coef_y = [+1, -2, +1]
	
	def __init__(self, pde = None):
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
		alpha = b_cond[0]
		beta  = b_cond[1]
		phi   = b_cond[2]
		N = len(u)
		h = self.h

		
		koef = (beta*h - alpha)
		
		return [0.0] + [(phi(i*h) * h - u[i]) / koef for i in range(1, N-1)] + [0.0]
	
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

	def Zeidel_method(self):
		"""Just solve equation"""
		U = deepcopy(self.grid[-1])
		N = self.N
		M = self.M
		h = self.h
		l = self.l
		
		tmp = [U[1][j] for j in range(M)] 
		top = boundary_eq(tmp, self.north, self.u_yy, self.u_xx, lambda z:f(z,0))
		tmp = [U[i][1] for i in range(N)]
		left = boundary_eq(tmp, self.west, self.u_xx, self.u_yy, lambda z:f(0,z))
		for j in range(1, M-1):
			U[0][j] = top[j]
		for i in range(1, N-1):
			U[i][0] = left[i]
		
		# U[i][j-1], U[i-1][j] takes from this iteration
		# U[i+1][j], U[i][j+1] takes from last iteration
		for j in range(1, M-1):
			for i in range(1, N-1):
				U1[i,j] = (U[i][j-1] + U[i-1][j] + U[i+1][j] + U[i][j+1] - h*h*f(i*h, j*h)) / 4
		
		tmp = [U1[-2][j] for j in range(M)] 
		footer = self.boundary_eq(tmp, self.south, self.u_yy, self.u_xx, lambda z:f(z,0))
		tmp = [U[i][-2] for i in range(N)]
		right = self.boundary_eq(tmp, self.east, self.u_xx, self.u_yy, lambda z:f(0,z))
		for j in range(1, M-1):
			U[-1][j] = footer[j]
		for i in range(1, N-1):
			U[i][-1] = right[i]

	def Liebmann_method(self, method = "Liebmann"):
		"""Just solve equation"""
		U = self.grid[-1]
		N = len(U)
		M = len(U[0])
		h = self.h
		lx = self.lx
		ly = self.ly
		fun = self.fun
		
		##correct size
		
		U = [[0.0]*N] + [[0.0] +
		      [(U[i][j-1] + U[i-1][j] + U[i+1][j] + U[i][j+1] - h*h * fun(i*h, j*h)) / 4
			  for j in range(1, M-1)] + [0.0]
			  for i in range(1, N-1)] + [[0.0]*N]
		
		#solve boundary conditions here!
		
		U[0] = self.boundary_eq(U[1], self.north, self.u_yy, self.u_xx, lambda z:fun(z,0))
		U[-1] = self.boundary_eq(U[-2], self.south, self.u_yy, self.u_xx, lambda z:fun(z,ly))
		U = [list(x) for x in zip(*U)] #transponate
		U[0] = self.boundary_eq(U[1], self.west, self.u_xx, self.u_yy, lambda z:fun(0,z))
		U[-1] = self.boundary_eq(U[-2], self.east, self.u_xx, self.u_yy, lambda z:fun(lx,z))
		U = [list(x) for x in zip(*U)] #transponate
		
		pprint(U)
		
		return U
		
	def estimate_error(self):
		U0 = self.grid[-1]
		U1 = self.grid[-2]
		
		eps = max(max(abs(e0-e1) for e0,e1 in zip(v0,v1))
		                         for v0,v1 in zip(U0,U1))
		return eps
	
	def solve(self, method = 'Liebmann'):
		"""Description"""
		self.M = int(self.lx / self.h)
		self.N = int(self.ly / self.h)
		
		self.grid = [[0.0]*self.M for i in range(self.N)]
		Us = [self.grid]
		self.grid = Us
		count = 0
		if method == 'Liebmann':
			Us.append(self.Liebmann_method())
			while self.estimate_error() > self.eps and count < 7:
				Us.append(self.Liebmann_method())
				if count == 2:
					print("RESULT:")
					pprint(Us[-1])
					return 
				count += 1
		elif method == 'Zeidel':
			while self.estimate_error() > self.eps and count < 7:
				Us.append(self.Zeidel_method())
				count += 1
		elif method == 'Relax':
			pass
		print("complete")
		return Us
	
	def check_scheme(self):
		""" GSOM GSOM """
		N = len(self.grid)
		M = len(self.grid[0])
		h = self.h
		fun = self.fun
		
		lx = self.lx
		ly = self.ly

		u = self.res_fun
		
		fu = [[u(i*h, j*h) for j in range(M)]
		                   for i in range(N)]

		errs = sum(sum(fu[i-1][j] - 2*fu[i][j] + fu[i+1][j] +
				fu[i][j-1] - 2*fu[i][j] + fu[i][j+1] -
				fun(i*h, j*h)
					for j in range(1,M-1))
					for i in range(1,N-1))
		
		
		alpha = self.north[0]
		beta = self.north[1]
		phi1 = self.north[2]
		
		print("let's have test u")
		print(u(0,0), u(0.5,0.5), u(0.5, 0.7))

		north_err = [alpha * (u(lx-h, h*j) - u(lx-0, h*j)) / h  +
		             beta *  u(lx-0, h*j) -
		             phi1(x = h*j) for j in range(M)]


		alpha = self.south[0]
		beta = self.south[1]
		phi0 = self.south[2]

		south_err = [alpha * (u(0, h*j)  - u(h, h*j)) / h  +
		             beta *  u(0, h*j) -
		             phi0(x = h*j) for j in range(M)]


		alpha = self.left[0]
		beta = self.left[1]
		phi0 = self.left[2]

		west_err = [alpha * (u(h*i, h) - u(h*i, 0)) / h  +
		             beta *  u(h*i, 0) -
		             phi0(y = h*i) for i in range(N)]

		alpha = self.right[0]
		beta = self.right[1]
		phi0 = self.right[2]

		east_err = [alpha * (u(h*i, ly-h) - u(h*i, ly-0)) / h  +
		            beta *  u(h*i, ly-0) -
		            phi0(y = h*i) for i in range(N)]


		
		return errs, north_err, south_err, west_err, east_err

	



#=====================================================================

def main():
	pass

if __name__ == "__main__":
	main()
