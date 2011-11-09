#! /usr/bin/env python3
# -*- coding: utf-8 -*-

from math import sqrt,sin,cos,tan,pi,log,exp,sqrt
from copy import copy, deepcopy
from functools import reduce
from tridiagonal import *

from pprint import pprint

from pde import *


class Hyperbolic_PDE(PDE):
	u_x = 0.0
	u = 0.0
	coef_a = [+1, -2, +1]
	coef_b = [-3, +4, -1]
	coef_c = [+0, +1, +0]
	coef_t = [+1, -2, +1]
	
	def __init__(self, pde = None):
		super(Parabolic_PDE, self).__init__(pde)
		MetaClass.print(self)
		
		a = sqrt(self.u_xx)
		b = self.u_x
		c = self.u
		h = self.h
		l = self.l
		tau = self.tau
		
		self.sigma = tau**2 * a**2 / h**2
		self.omega = tau**2 * b / (2*h)
		self.eta = tau**2 * c # TODO add f(x,t)

		self.coefficients = PDE.vec_mat([self.sigma, self.omega, self.eta],
	                       [self.coef_a, self.coef_b, self.coef_c])
	    
		psi0 = self.initial0
		U = []
		psi0 = self.initial0
		psi1 = self.initial1
		
		U.append([ psi0(x) for x in frange(0, T, t)])
		U.append([ psi0(x) + psi1(x)*tau for x in frange(0, T, t)])
		for k in range(1, N-1):
			U[-1][k] += (U[0][k-1] - 2*U[0][k] + U[0][k+1]) * alpha**2 * tau**2 / (2 * h**2)
 		
		self.grid = U
		

	

	def middle_eq(self, i, teta = 1):
		"""Find coefficients of middle equation"""
		U = self.grid[-1]
		N = len(U)
		
		coeff = PDE.vec_mat([self.sigma, self.omega, self.eta],
		                    [self.coef_a, self.coef_b, self.coef_c])
	
		
		ai = teta * coeff[0]
		bi = teta * coeff[1] - 1
		ci = teta * coeff[2]
		di = -U[i]
		
		
		return (ai, bi, ci, di)
		
	def last_eq(u):
		"""Find coefficients of last equation"""
		pass

		
	def solve(self):
		U = self.grid
		print("ololo")

#=====================================================================

def main():

	f = open("input")
	
	pde = parse_file(f)

#	print(pde.left[2](0))
#	print(pde.right[2](0))
	
	res = pde.solve()
	
	for k in [0,1,2,3]:
		print("iter ", k)
		print_vec(res[k])
		print_vec([pde.fun(x, k*pde.tau) for x in frange(0, pde.l, pde.h)])
	
#	print(pde.u_x)


#=====================================================================

if __name__ == "__main__":
	main()
